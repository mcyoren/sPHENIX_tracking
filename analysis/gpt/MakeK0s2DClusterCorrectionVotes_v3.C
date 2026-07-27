// MakeK0s2DClusterCorrectionVotes_v3.C
//
// K0S -> pi+ pi- constrained transverse cluster-correction votes.
//
// For every clean low-pT K0S candidate:
//   1. keep each daughter trajectory anchored near R2;
//   2. alternate momentum-scale and small-rotation fits;
//   3. constrain m(pi pi) to the PDG K0S mass, DIRA toward 1, and pair DCA to 0;
//   4. compare each measured cluster with the corrected daughter circle;
//   5. fill (Delta r, r Delta phi) vote histograms.
//
// Final detector layout:
//   modules  0-11 : R1 sectors 0-11
//   modules 12-23 : R2 sectors 0-11
//   modules 24-35 : R3 sectors 0-11
//
// Vote directory:
//   side{0,1}/{qplus,qminus}/pt_.../module_XX/phi_I_layer_J/
//       h_deltaRPhi_vs_deltaR_votes
//
// The companion notebook normalizes each charge-pT category before combining
// them, so the final correction map is pT weighted rather than dominated by
// the highest-statistics category.
//
// Example:
// root -l -b -q 'MakeK0s2DClusterCorrectionVotes_v3.C(
//   "/path/to/pair/files","p_v0_*.root",
//   "k0s_2d_cluster_votes.root","pairTree")'

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
  constexpr double kPionMass = 0.13957039;
  constexpr double kKshortMass = 0.497611;
  constexpr double kTwoPi = 2.0 * TMath::Pi();

  struct Vec3
  {
    double x = 0.;
    double y = 0.;
    double z = 0.;
  };

  struct Circle
  {
    double xc = 0.;
    double yc = 0.;
    double radius = 0.;
    bool valid = false;
  };

  struct Daughter
  {
    Int_t charge = 0;
    Int_t side = -1;
    Int_t npoints = 0;
    UInt_t ntpcClusters = 0;

    Float_t px = 0.F;
    Float_t py = 0.F;
    Float_t pz = 0.F;
    Float_t pt = 0.F;
    Float_t eta = 0.F;

    std::vector<unsigned int>* layer = nullptr;
    std::vector<double>* clusterX = nullptr;
    std::vector<double>* clusterY = nullptr;
    std::vector<double>* clusterZ = nullptr;
    std::vector<double>* clusterR = nullptr;
    std::vector<double>* clusterPhi = nullptr;

    std::vector<double>* fitX = nullptr;
    std::vector<double>* fitY = nullptr;
    std::vector<double>* fitZ = nullptr;
    std::vector<double>* fitPx = nullptr;
    std::vector<double>* fitPy = nullptr;
    std::vector<double>* fitPz = nullptr;
  };

  struct Anchor
  {
    Vec3 position;
    Vec3 momentum;
    bool valid = false;
  };

  struct FitResult
  {
    bool valid = false;
    int iterations = 0;

    double scale1 = 1.;
    double scale2 = 1.;
    double rotation1 = 0.;
    double rotation2 = 0.;

    double massBefore = 0.;
    double massAfter = 0.;

    double diraXYBefore = -2.;
    double diraXYAfter = -2.;
    double dira3DBefore = -2.;
    double dira3DAfter = -2.;

    double pairDcaXYBefore = 0.;
    double pairDcaXYAfter = 0.;
    double pairDca3DBefore = 0.;
    double pairDca3DAfter = 0.;

    Vec3 momentum1;
    Vec3 momentum2;
    Circle circle1;
    Circle circle2;
  };


  struct V0Cut
  {
    std::string name;
    double maximumAbsPcaZ = 0.;
    double maximumAbsDeltaPcaZ = 0.;
    double minimumDaughterPt = 0.;
    double minimumDecayRadius = 0.;
    double maximumAbsAlpha = 0.;
    double maximumPairDca = 0.;
    double minimumDira = 0.;
    double maximumQuality = 0.;
    int minimumNpoints = 0;
  };

  V0Cut makeRequestedCut(const std::string& name)
  {
    if (name == "cut03_baseline")
    {
      return {"cut03_baseline",15.0,0.50,0.20,2.0,0.99,2.00,0.85,15.0,30};
    }

    if (name == "cut07_pairDCA_5mm")
    {
      return {"cut07_pairDCA_5mm",10.0,0.20,0.20,2.0,0.99,0.50,0.95,10.0,32};
    }

    std::cerr << "WARNING: unknown selectionName '" << name
              << "'. Using cut03_baseline." << std::endl;
    return makeRequestedCut("cut03_baseline");
  }

  Vec3 add(const Vec3& a, const Vec3& b)
  {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  }

  Vec3 subtract(const Vec3& a, const Vec3& b)
  {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  }

  Vec3 multiply(const Vec3& a, const double value)
  {
    return {value * a.x, value * a.y, value * a.z};
  }

  double dot(const Vec3& a, const Vec3& b)
  {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  double norm2(const Vec3& a)
  {
    return dot(a, a);
  }

  double norm(const Vec3& a)
  {
    return std::sqrt(norm2(a));
  }

  double transverseMomentum(const Vec3& a)
  {
    return std::hypot(a.x, a.y);
  }

  bool finite(const Vec3& a)
  {
    return std::isfinite(a.x) &&
           std::isfinite(a.y) &&
           std::isfinite(a.z);
  }

  Vec3 unit(const Vec3& a)
  {
    const double magnitude = norm(a);
    if (!(magnitude > 0.) || !std::isfinite(magnitude))
    {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      return {nan, nan, nan};
    }
    return multiply(a, 1. / magnitude);
  }

  Vec3 rotateZ(const Vec3& p, const double angle)
  {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return {
      c * p.x - s * p.y,
      s * p.x + c * p.y,
      p.z
    };
  }

  double wrapPhi(double phi)
  {
    while (phi >= TMath::Pi()) phi -= kTwoPi;
    while (phi < -TMath::Pi()) phi += kTwoPi;
    return phi;
  }

  double wrapToTwoPi(double phi)
  {
    phi = std::fmod(phi, kTwoPi);
    if (phi < 0.) phi += kTwoPi;
    return phi;
  }

  double invariantMass(const Vec3& p1, const Vec3& p2)
  {
    const double e1 = std::sqrt(norm2(p1) + kPionMass * kPionMass);
    const double e2 = std::sqrt(norm2(p2) + kPionMass * kPionMass);
    const Vec3 total = add(p1, p2);
    const double mass2 = (e1 + e2) * (e1 + e2) - norm2(total);
    return mass2 > 0. ? std::sqrt(mass2) : 0.;
  }

  double calculateDira(const Vec3& momentum, const Vec3& flight)
  {
    const double denominator = norm(momentum) * norm(flight);
    if (!(denominator > 0.)) return -2.;
    return std::clamp(dot(momentum, flight) / denominator, -1., 1.);
  }

  double openingAngle(const Vec3& p1, const Vec3& p2)
  {
    const double denominator = norm(p1) * norm(p2);
    if (!(denominator > 0.))
      return std::numeric_limits<double>::quiet_NaN();

    const double cosine =
        std::clamp(dot(p1, p2) / denominator, -1., 1.);
    return std::acos(cosine);
  }

  double lineLineDca(
      const Vec3& point1,
      const Vec3& direction1Input,
      const Vec3& point2,
      const Vec3& direction2Input)
  {
    const Vec3 direction1 = unit(direction1Input);
    const Vec3 direction2 = unit(direction2Input);

    if (!finite(direction1) || !finite(direction2))
      return std::numeric_limits<double>::quiet_NaN();

    const Vec3 difference = subtract(point1, point2);

    const double a = dot(direction1, direction1);
    const double b = dot(direction1, direction2);
    const double c = dot(direction2, direction2);
    const double d = dot(direction1, difference);
    const double e = dot(direction2, difference);
    const double denominator = a * c - b * b;

    double parameter1 = 0.;
    double parameter2 = 0.;

    if (std::abs(denominator) > 1e-12)
    {
      parameter1 = (b * e - c * d) / denominator;
      parameter2 = (a * e - b * d) / denominator;
    }
    else
    {
      parameter2 = c > 0. ? e / c : 0.;
    }

    const Vec3 closest1 =
        add(point1, multiply(direction1, parameter1));
    const Vec3 closest2 =
        add(point2, multiply(direction2, parameter2));

    return norm(subtract(closest1, closest2));
  }

  int ptBin(const double pt)
  {
    if (pt >= 0.2 && pt < 0.4) return 0;
    if (pt >= 0.4 && pt < 0.7) return 1;
    if (pt >= 0.7 && pt < 1.2) return 2;
    if (pt >= 1.2 && pt < 1.8) return 3;
    if (pt >= 1.8 && pt < 5.0) return 4;
    return -1;
  }

  std::string ptBinName(const int bin)
  {
    static const std::array<const char*, 5> names = {
      "pt_0p2_0p4",
      "pt_0p4_0p7",
      "pt_0p7_1p2",
      "pt_1p2_1p8",
      "pt_1p8_5p0"
    };

    return (bin >= 0 && bin < static_cast<int>(names.size()))
      ? names[bin]
      : "pt_invalid";
  }

  int radialModuleFromLayer(const int layer)
  {
    if (layer >= 7 && layer <= 22) return 0;
    if (layer >= 23 && layer <= 38) return 1;
    if (layer >= 39 && layer <= 54) return 2;
    return -1;
  }

  int localLayerBin(const int layer, const int radialModule)
  {
    if (radialModule < 0 || radialModule > 2) return -1;

    const int firstLayer = 7 + 16 * radialModule;
    const int localLayer = layer - firstLayer;

    // Remove first and last layer of each radial region.
    if (localLayer <= 0 || localLayer >= 15) return -1;

    const int index = localLayer - 1;
    if (index < 3) return 0;
    if (index < 6) return 1;
    if (index < 9) return 2;
    if (index < 12) return 3;
    return 4;
  }

  struct PhiBin
  {
    int sector = -1;
    int localPhi = -1;
  };

  PhiBin determinePhiBin(
      const double phi,
      const double sectorPhiOffset)
  {
    PhiBin result;

    const double sectorWidth = kTwoPi / 12.;
    const double shifted = wrapToTwoPi(phi - sectorPhiOffset);

    result.sector =
        std::clamp(static_cast<int>(shifted / sectorWidth), 0, 11);

    const double withinSector =
        shifted - result.sector * sectorWidth;
    const double localFraction = withinSector / sectorWidth;

    result.localPhi =
        std::clamp(static_cast<int>(3. * localFraction), 0, 2);

    return result;
  }

  int hardwareModule(
      const int sector,
      const int radialModule)
  {
    // Explicit requested convention:
    // R1: 0-11, R2: 12-23, R3: 24-35.
    if (sector < 0 || sector > 11 ||
        radialModule < 0 || radialModule > 2)
    {
      return -1;
    }

    return 12 * radialModule + sector;
  }

  std::string histogramPath(
      const std::string& selectionName,
      const int side,
      const int charge,
      const int iPt,
      const int module,
      const int localPhi,
      const int localLayer)
  {
    return TString::Format(
      "%s/side%d/%s/%s/module_%02d/phi_%d_layer_%d",
      selectionName.c_str(),
      side,
      charge > 0 ? "qplus" : "qminus",
      ptBinName(iPt).c_str(),
      module,
      localPhi,
      localLayer).Data();
  }

  TDirectory* makeDirectoryPath(
      TDirectory* root,
      const std::string& path)
  {
    TDirectory* current = root;
    TString temporary(path.c_str());
    std::unique_ptr<TObjArray> tokens(temporary.Tokenize("/"));

    for (int index = 0; index < tokens->GetEntries(); ++index)
    {
      const TString token = tokens->At(index)->GetName();
      TDirectory* next = current->GetDirectory(token);
      if (!next) next = current->mkdir(token);
      current = next;
    }

    return current;
  }

  TH2D* bookVoteHistogram(
      TDirectory* directory,
      const std::string& label,
      const int bins,
      const double maximumCorrection)
  {
    directory->cd();

    TH2D* histogram = new TH2D(
      "h_deltaRPhi_vs_deltaR_votes",
      TString::Format(
        "K^{0}_{S}-constrained cluster votes [%s];"
        "#Delta r [cm];r#Delta#phi [cm]",
        label.c_str()),
      bins,
      -maximumCorrection,
      maximumCorrection,
      bins,
      -maximumCorrection,
      maximumCorrection);

    histogram->Sumw2();
    return histogram;
  }

  void bindDaughter(
      TChain& chain,
      const std::string& prefix,
      Daughter& daughter)
  {
    chain.SetBranchAddress(
      (prefix + "_charge").c_str(),
      &daughter.charge);
    chain.SetBranchAddress(
      (prefix + "_side").c_str(),
      &daughter.side);
    chain.SetBranchAddress(
      (prefix + "_npoints").c_str(),
      &daughter.npoints);
    chain.SetBranchAddress(
      (prefix + "_ntpc_clusters").c_str(),
      &daughter.ntpcClusters);

    chain.SetBranchAddress(
      (prefix + "_px").c_str(),
      &daughter.px);
    chain.SetBranchAddress(
      (prefix + "_py").c_str(),
      &daughter.py);
    chain.SetBranchAddress(
      (prefix + "_pz").c_str(),
      &daughter.pz);
    chain.SetBranchAddress(
      (prefix + "_pt").c_str(),
      &daughter.pt);
    chain.SetBranchAddress(
      (prefix + "_eta").c_str(),
      &daughter.eta);

    chain.SetBranchAddress(
      (prefix + "_layer").c_str(),
      &daughter.layer);
    chain.SetBranchAddress(
      (prefix + "_cluster_x").c_str(),
      &daughter.clusterX);
    chain.SetBranchAddress(
      (prefix + "_cluster_y").c_str(),
      &daughter.clusterY);
    chain.SetBranchAddress(
      (prefix + "_cluster_z").c_str(),
      &daughter.clusterZ);
    chain.SetBranchAddress(
      (prefix + "_cluster_r").c_str(),
      &daughter.clusterR);
    chain.SetBranchAddress(
      (prefix + "_cluster_phi").c_str(),
      &daughter.clusterPhi);

    chain.SetBranchAddress(
      (prefix + "_fit_x").c_str(),
      &daughter.fitX);
    chain.SetBranchAddress(
      (prefix + "_fit_y").c_str(),
      &daughter.fitY);
    chain.SetBranchAddress(
      (prefix + "_fit_z").c_str(),
      &daughter.fitZ);
    chain.SetBranchAddress(
      (prefix + "_fit_px").c_str(),
      &daughter.fitPx);
    chain.SetBranchAddress(
      (prefix + "_fit_py").c_str(),
      &daughter.fitPy);
    chain.SetBranchAddress(
      (prefix + "_fit_pz").c_str(),
      &daughter.fitPz);
  }

  bool daughterVectorsValid(const Daughter& daughter)
  {
    return daughter.layer &&
           daughter.clusterX &&
           daughter.clusterY &&
           daughter.clusterZ &&
           daughter.clusterR &&
           daughter.clusterPhi &&
           daughter.fitX &&
           daughter.fitY &&
           daughter.fitZ &&
           daughter.fitPx &&
           daughter.fitPy &&
           daughter.fitPz;
  }

  std::size_t daughterVectorSize(const Daughter& daughter)
  {
    if (!daughterVectorsValid(daughter)) return 0;

    return std::min({
      daughter.layer->size(),
      daughter.clusterX->size(),
      daughter.clusterY->size(),
      daughter.clusterZ->size(),
      daughter.clusterR->size(),
      daughter.clusterPhi->size(),
      daughter.fitX->size(),
      daughter.fitY->size(),
      daughter.fitZ->size(),
      daughter.fitPx->size(),
      daughter.fitPy->size(),
      daughter.fitPz->size()
    });
  }

  Anchor chooseR2Anchor(
      const Daughter& daughter,
      const double anchorMinimum,
      const double anchorMaximum,
      const double preferredRadius)
  {
    Anchor result;

    const std::size_t size = daughterVectorSize(daughter);
    double bestDistance = std::numeric_limits<double>::max();

    for (std::size_t index = 0; index < size; ++index)
    {
      const double radius = daughter.clusterR->at(index);
      if (!std::isfinite(radius) ||
          radius < anchorMinimum ||
          radius > anchorMaximum)
      {
        continue;
      }

      const Vec3 position{
        daughter.fitX->at(index),
        daughter.fitY->at(index),
        daughter.fitZ->at(index)
      };

      const Vec3 momentum{
        daughter.fitPx->at(index),
        daughter.fitPy->at(index),
        daughter.fitPz->at(index)
      };

      if (!finite(position) || !finite(momentum) ||
          !(transverseMomentum(momentum) > 0.))
      {
        continue;
      }

      const double distance =
          std::abs(radius - preferredRadius);

      if (distance < bestDistance)
      {
        bestDistance = distance;
        result.position = position;
        result.momentum = momentum;
        result.valid = true;
      }
    }

    return result;
  }

  Circle circleFromAnchor(
      const Anchor& anchor,
      const double momentumScale,
      const double rotation,
      const int charge,
      const double bFieldTesla)
  {
    Circle result;
    if (!anchor.valid || std::abs(charge) < 1) return result;

    const Vec3 corrected =
        rotateZ(multiply(anchor.momentum, momentumScale), rotation);

    const double correctedPt = transverseMomentum(corrected);
    const double denominator =
        0.003 * static_cast<double>(charge) * bFieldTesla;

    if (!(correctedPt > 0.) ||
        !std::isfinite(denominator) ||
        std::abs(denominator) < 1e-12)
    {
      return result;
    }

    const double signedRadius = correctedPt / denominator;
    const double tangentX = corrected.x / correctedPt;
    const double tangentY = corrected.y / correctedPt;

    result.xc =
        anchor.position.x + signedRadius * tangentY;
    result.yc =
        anchor.position.y - signedRadius * tangentX;
    result.radius = std::abs(signedRadius);
    result.valid =
        std::isfinite(result.xc) &&
        std::isfinite(result.yc) &&
        std::isfinite(result.radius) &&
        result.radius > 0.;

    return result;
  }

  std::pair<double, double> closestPointOnCircle(
      const Circle& circle,
      const double x,
      const double y)
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!circle.valid) return {nan, nan};

    const double dx = x - circle.xc;
    const double dy = y - circle.yc;
    const double distance = std::hypot(dx, dy);

    if (!(distance > 0.) || !std::isfinite(distance))
      return {nan, nan};

    return {
      circle.xc + circle.radius * dx / distance,
      circle.yc + circle.radius * dy / distance
    };
  }

  double calculateDiraXY(
      const Vec3& momentum,
      const Vec3& flight)
  {
    const double momentumMagnitude =
        std::hypot(momentum.x, momentum.y);
    const double flightMagnitude =
        std::hypot(flight.x, flight.y);

    if (!(momentumMagnitude > 0.) ||
        !(flightMagnitude > 0.))
    {
      return -2.;
    }

    return std::clamp(
      (momentum.x * flight.x +
       momentum.y * flight.y) /
      (momentumMagnitude * flightMagnitude),
      -1.,
      1.);
  }

  double transverseCross(
      const Vec3& first,
      const Vec3& second)
  {
    return first.x * second.y -
           first.y * second.x;
  }

  double circleCircleDcaXY(
      const Circle& first,
      const Circle& second)
  {
    if (!first.valid || !second.valid)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    const double centerDistance =
        std::hypot(
          first.xc - second.xc,
          first.yc - second.yc);

    if (!std::isfinite(centerDistance))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    if (centerDistance > first.radius + second.radius)
    {
      return centerDistance -
             first.radius -
             second.radius;
    }

    if (centerDistance < std::abs(first.radius - second.radius))
    {
      return std::abs(first.radius - second.radius) -
             centerDistance;
    }

    // The circles intersect in XY.
    return 0.;
  }

  struct ScaleSolution
  {
    bool valid = false;
    double scale1 = 1.;
    double scale2 = 1.;
    Vec3 momentum1;
    Vec3 momentum2;
    double mass = 0.;
    double diraXY = -2.;
  };

  ScaleSolution solveExactTransverseScales(
      const Vec3& originalMomentum1,
      const Vec3& originalMomentum2,
      const double rotation1,
      const double rotation2,
      const Vec3& flight,
      const double minimumScale,
      const double maximumScale)
  {
    ScaleSolution result;

    const Vec3 rotated1 =
        rotateZ(originalMomentum1, rotation1);
    const Vec3 rotated2 =
        rotateZ(originalMomentum2, rotation2);

    const double denominator =
        transverseCross(rotated2, flight);

    if (!finite(rotated1) ||
        !finite(rotated2) ||
        std::abs(denominator) < 1e-12)
    {
      return result;
    }

    // Exact transverse pointing:
    // s1 * cross(p1,L) + s2 * cross(p2,L) = 0.
    const double ratio =
        -transverseCross(rotated1, flight) /
         denominator;

    if (!(ratio > 0.) || !std::isfinite(ratio))
    {
      return result;
    }

    const double lower =
        std::max(minimumScale, minimumScale / ratio);
    const double upper =
        std::min(maximumScale, maximumScale / ratio);

    if (!(upper > lower))
    {
      return result;
    }

    auto massAt = [&](const double commonScale)
    {
      const Vec3 momentum1 =
          multiply(rotated1, commonScale);
      const Vec3 momentum2 =
          multiply(rotated2, ratio * commonScale);
      return invariantMass(momentum1, momentum2);
    };

    // Find a bracket for m(commonScale) - mPDG.
    constexpr int scanSteps = 500;
    double previousScale = lower;
    double previousValue = massAt(previousScale) - kKshortMass;

    bool bracketFound = false;
    double bracketLow = lower;
    double bracketHigh = upper;

    double bestScale = lower;
    double bestResidual = std::abs(previousValue);

    for (int index = 1; index <= scanSteps; ++index)
    {
      const double fraction =
          static_cast<double>(index) / scanSteps;
      const double currentScale =
          lower + fraction * (upper - lower);
      const double currentValue =
          massAt(currentScale) - kKshortMass;

      if (std::abs(currentValue) < bestResidual)
      {
        bestResidual = std::abs(currentValue);
        bestScale = currentScale;
      }

      if (std::isfinite(previousValue) &&
          std::isfinite(currentValue) &&
          previousValue * currentValue <= 0.)
      {
        bracketLow = previousScale;
        bracketHigh = currentScale;
        bracketFound = true;
        break;
      }

      previousScale = currentScale;
      previousValue = currentValue;
    }

    double commonScale = bestScale;

    if (bracketFound)
    {
      double low = bracketLow;
      double high = bracketHigh;
      double valueLow = massAt(low) - kKshortMass;

      for (int iteration = 0; iteration < 80; ++iteration)
      {
        const double middle = 0.5 * (low + high);
        const double valueMiddle =
            massAt(middle) - kKshortMass;

        if (!std::isfinite(valueMiddle))
        {
          return result;
        }

        commonScale = middle;

        if (std::abs(valueMiddle) < 1e-10)
        {
          break;
        }

        if (valueLow * valueMiddle <= 0.)
        {
          high = middle;
        }
        else
        {
          low = middle;
          valueLow = valueMiddle;
        }
      }
    }

    result.scale1 = commonScale;
    result.scale2 = ratio * commonScale;
    result.momentum1 =
        multiply(rotated1, result.scale1);
    result.momentum2 =
        multiply(rotated2, result.scale2);
    result.mass =
        invariantMass(result.momentum1, result.momentum2);
    result.diraXY =
        calculateDiraXY(
          add(result.momentum1, result.momentum2),
          flight);

    const Vec3 total =
        add(result.momentum1, result.momentum2);
    const double forward =
        total.x * flight.x +
        total.y * flight.y;

    result.valid =
        finite(result.momentum1) &&
        finite(result.momentum2) &&
        std::isfinite(result.mass) &&
        std::isfinite(result.diraXY) &&
        forward > 0. &&
        std::abs(result.mass - kKshortMass) < 1e-5 &&
        result.diraXY > 1. - 1e-8 &&
        result.scale1 >= minimumScale &&
        result.scale1 <= maximumScale &&
        result.scale2 >= minimumScale &&
        result.scale2 <= maximumScale;

    return result;
  }

  struct RotationSolution
  {
    bool valid = false;
    double rotation1 = 0.;
    double rotation2 = 0.;
    double pairDcaXY = 0.;
    Circle circle1;
    Circle circle2;
  };

  RotationSolution solveRotationsFromPairDcaXY(
      const Anchor& anchor1,
      const Anchor& anchor2,
      const int charge1,
      const int charge2,
      const double scale1,
      const double scale2,
      const double magneticFieldTesla,
      const double currentRotation1,
      const double currentRotation2,
      const double maximumRotation)
  {
    RotationSolution result;

    double center1 = currentRotation1;
    double center2 = currentRotation2;
    double halfWidth = maximumRotation;

    for (int refinement = 0; refinement < 5; ++refinement)
    {
      constexpr int steps = 60;

      double bestObjective =
          std::numeric_limits<double>::max();
      double bestRotation1 = center1;
      double bestRotation2 = center2;
      double bestDca =
          std::numeric_limits<double>::quiet_NaN();
      Circle bestCircle1;
      Circle bestCircle2;

      for (int index1 = 0; index1 <= steps; ++index1)
      {
        const double fraction1 =
            static_cast<double>(index1) / steps;

        const double candidateRotation1 =
            std::clamp(
              center1 - halfWidth +
                2. * halfWidth * fraction1,
              -maximumRotation,
              maximumRotation);

        for (int index2 = 0; index2 <= steps; ++index2)
        {
          const double fraction2 =
              static_cast<double>(index2) / steps;

          const double candidateRotation2 =
              std::clamp(
                center2 - halfWidth +
                  2. * halfWidth * fraction2,
                -maximumRotation,
                maximumRotation);

          const Circle circle1 =
              circleFromAnchor(
                anchor1,
                scale1,
                candidateRotation1,
                charge1,
                magneticFieldTesla);

          const Circle circle2 =
              circleFromAnchor(
                anchor2,
                scale2,
                candidateRotation2,
                charge2,
                magneticFieldTesla);

          const double dca =
              circleCircleDcaXY(circle1, circle2);

          if (!std::isfinite(dca)) continue;

          // pairDCA_xy is the primary target.  The very small
          // rotation term only selects the least-rotated solution
          // when many circle pairs intersect.
          const double rotationNorm =
              candidateRotation1 * candidateRotation1 +
              candidateRotation2 * candidateRotation2;

          const double objective =
              dca * dca +
              1e-6 * rotationNorm;

          if (objective < bestObjective)
          {
            bestObjective = objective;
            bestRotation1 = candidateRotation1;
            bestRotation2 = candidateRotation2;
            bestDca = dca;
            bestCircle1 = circle1;
            bestCircle2 = circle2;
          }
        }
      }

      if (!std::isfinite(bestDca))
      {
        return result;
      }

      center1 = bestRotation1;
      center2 = bestRotation2;
      result.rotation1 = bestRotation1;
      result.rotation2 = bestRotation2;
      result.pairDcaXY = bestDca;
      result.circle1 = bestCircle1;
      result.circle2 = bestCircle2;
      halfWidth *= 0.20;
    }

    result.valid =
        result.circle1.valid &&
        result.circle2.valid &&
        std::isfinite(result.pairDcaXY);

    return result;
  }

  FitResult solveAlternating(
      const Vec3& originalMomentum1,
      const Vec3& originalMomentum2,
      const Vec3& flight,
      const Vec3& pca1,
      const Vec3& pca2,
      const Anchor& anchor1,
      const Anchor& anchor2,
      const int charge1,
      const int charge2,
      const double magneticFieldTesla,
      const int requestedIterations,
      const double minimumScale,
      const double maximumScale,
      const double maximumRotation)
  {
    FitResult result;

    result.massBefore =
        invariantMass(originalMomentum1, originalMomentum2);
    result.diraXYBefore =
        calculateDiraXY(
          add(originalMomentum1, originalMomentum2),
          flight);
    result.dira3DBefore =
        calculateDira(
          add(originalMomentum1, originalMomentum2),
          flight);
    result.pairDca3DBefore =
        lineLineDca(
          pca1, originalMomentum1,
          pca2, originalMomentum2);

    const Circle originalCircle1 =
        circleFromAnchor(
          anchor1,
          1.0,
          0.0,
          charge1,
          magneticFieldTesla);

    const Circle originalCircle2 =
        circleFromAnchor(
          anchor2,
          1.0,
          0.0,
          charge2,
          magneticFieldTesla);

    result.pairDcaXYBefore =
        circleCircleDcaXY(
          originalCircle1,
          originalCircle2);

    double rotation1 = 0.;
    double rotation2 = 0.;
    ScaleSolution scales;
    RotationSolution rotations;

    for (int iteration = 0;
         iteration < requestedIterations;
         ++iteration)
    {
      scales = solveExactTransverseScales(
        originalMomentum1,
        originalMomentum2,
        rotation1,
        rotation2,
        flight,
        minimumScale,
        maximumScale);

      if (!scales.valid)
      {
        return result;
      }

      rotations = solveRotationsFromPairDcaXY(
        anchor1,
        anchor2,
        charge1,
        charge2,
        scales.scale1,
        scales.scale2,
        magneticFieldTesla,
        rotation1,
        rotation2,
        maximumRotation);

      if (!rotations.valid)
      {
        return result;
      }

      rotation1 = rotations.rotation1;
      rotation2 = rotations.rotation2;
      result.iterations = iteration + 1;
    }

    // Final exact mass + transverse-pointing solution after
    // the last rotation update.
    scales = solveExactTransverseScales(
      originalMomentum1,
      originalMomentum2,
      rotation1,
      rotation2,
      flight,
      minimumScale,
      maximumScale);

    if (!scales.valid)
    {
      return result;
    }

    result.scale1 = scales.scale1;
    result.scale2 = scales.scale2;
    result.rotation1 = rotation1;
    result.rotation2 = rotation2;
    result.momentum1 = scales.momentum1;
    result.momentum2 = scales.momentum2;

    result.circle1 =
        circleFromAnchor(
          anchor1,
          result.scale1,
          result.rotation1,
          charge1,
          magneticFieldTesla);

    result.circle2 =
        circleFromAnchor(
          anchor2,
          result.scale2,
          result.rotation2,
          charge2,
          magneticFieldTesla);

    result.massAfter =
        invariantMass(result.momentum1, result.momentum2);
    result.diraXYAfter =
        calculateDiraXY(
          add(result.momentum1, result.momentum2),
          flight);
    result.dira3DAfter =
        calculateDira(
          add(result.momentum1, result.momentum2),
          flight);
    result.pairDcaXYAfter =
        circleCircleDcaXY(
          result.circle1,
          result.circle2);
    result.pairDca3DAfter =
        lineLineDca(
          pca1, result.momentum1,
          pca2, result.momentum2);

    result.valid =
        result.circle1.valid &&
        result.circle2.valid &&
        std::isfinite(result.massAfter) &&
        std::isfinite(result.diraXYAfter) &&
        std::isfinite(result.dira3DAfter) &&
        std::isfinite(result.pairDcaXYAfter) &&
        std::isfinite(result.pairDca3DAfter) &&
        std::abs(result.massAfter - kKshortMass) < 1e-5 &&
        result.diraXYAfter > 1. - 1e-8;

    return result;
  }

  // Fill the complete visible correction line
  //
  //   n_r Delta r + n_phi (r Delta phi) = d_normal
  //
  // with the requested total weight.
  void fillAllowedCorrectionLine(
      TH2D* histogram,
      const double normalR,
      const double normalPhi,
      const double normalDistance,
      const double maximumCorrection,
      const int lineSamples,
      const double totalWeight)
  {
    if (!histogram ||
        lineSamples < 2 ||
        !(totalWeight > 0.))
    {
      return;
    }

    std::vector<std::pair<double, double>> points;
    points.reserve(lineSamples);

    if (std::abs(normalPhi) >= std::abs(normalR) &&
        std::abs(normalPhi) > 1e-10)
    {
      for (int index = 0; index < lineSamples; ++index)
      {
        const double deltaR =
            -maximumCorrection +
            2. * maximumCorrection *
              static_cast<double>(index) /
              (lineSamples - 1);

        const double deltaRPhi =
            (normalDistance -
             normalR * deltaR) /
            normalPhi;

        if (std::abs(deltaRPhi) <= maximumCorrection)
        {
          points.emplace_back(deltaR, deltaRPhi);
        }
      }
    }
    else if (std::abs(normalR) > 1e-10)
    {
      for (int index = 0; index < lineSamples; ++index)
      {
        const double deltaRPhi =
            -maximumCorrection +
            2. * maximumCorrection *
              static_cast<double>(index) /
              (lineSamples - 1);

        const double deltaR =
            (normalDistance -
             normalPhi * deltaRPhi) /
            normalR;

        if (std::abs(deltaR) <= maximumCorrection)
        {
          points.emplace_back(deltaR, deltaRPhi);
        }
      }
    }

    if (points.empty()) return;

    const double pointWeight =
        totalWeight /
        static_cast<double>(points.size());

    for (const auto& point : points)
    {
      histogram->Fill(
        point.first,
        point.second,
        pointWeight);
    }
  }

  std::size_t countFillableClusters(
      const Daughter& daughter,
      const double maximumAbsClusterZ,
      const double maximumCorrection,
      const Circle& circle,
      const double sectorPhiOffset)
  {
    std::size_t count = 0;
    const std::size_t size = daughterVectorSize(daughter);

    for (std::size_t index = 0; index < size; ++index)
    {
      const int layer =
          static_cast<int>(daughter.layer->at(index));
      const int radialModule =
          radialModuleFromLayer(layer);
      const int localLayer =
          localLayerBin(layer, radialModule);

      if (radialModule < 0 || localLayer < 0) continue;

      const double x = daughter.clusterX->at(index);
      const double y = daughter.clusterY->at(index);
      const double z = daughter.clusterZ->at(index);
      const double radius = daughter.clusterR->at(index);
      const double phi = wrapPhi(daughter.clusterPhi->at(index));

      if (!std::isfinite(x) ||
          !std::isfinite(y) ||
          !std::isfinite(z) ||
          !std::isfinite(radius) ||
          !std::isfinite(phi) ||
          radius <= 0. ||
          std::abs(z) >= maximumAbsClusterZ)
      {
        continue;
      }

      const PhiBin phiBin =
          determinePhiBin(phi, sectorPhiOffset);
      const int module =
          hardwareModule(phiBin.sector, radialModule);

      if (module < 0 || phiBin.localPhi < 0) continue;

      const auto corrected =
          closestPointOnCircle(circle, x, y);

      if (!std::isfinite(corrected.first) ||
          !std::isfinite(corrected.second))
      {
        continue;
      }

      const double correctionX = corrected.first - x;
      const double correctionY = corrected.second - y;

      const double radialX = std::cos(phi);
      const double radialY = std::sin(phi);
      const double tangentialX = -std::sin(phi);
      const double tangentialY = std::cos(phi);

      const double deltaR =
          correctionX * radialX +
          correctionY * radialY;
      const double deltaRPhi =
          correctionX * tangentialX +
          correctionY * tangentialY;

      if (!std::isfinite(deltaR) ||
          !std::isfinite(deltaRPhi) ||
          std::abs(deltaR) > maximumCorrection ||
          std::abs(deltaRPhi) > maximumCorrection)
      {
        continue;
      }

      ++count;
    }

    return count;
  }

  void fillDaughterVotes(
      const std::string& selectionName,
      const Daughter& daughter,
      const double originalPt,
      const Circle& correctedCircle,
      TFile* output,
      std::map<std::string, TH2D*>& voteHistograms,
      std::map<std::string, Long64_t>& voteCounts,
      const double maximumAbsClusterZ,
      const double maximumCorrection,
      const int voteBins,
      const int lineSamples,
      const double sectorPhiOffset,
      const bool normalizeClustersPerDaughter,
      Long64_t& selectedClusters)
  {
    const int iPt = ptBin(originalPt);
    if (iPt < 0 || !correctedCircle.valid) return;

    const std::size_t fillableClusters =
        countFillableClusters(
          daughter,
          maximumAbsClusterZ,
          maximumCorrection,
          correctedCircle,
          sectorPhiOffset);

    if (fillableClusters == 0) return;

    const double weight =
        normalizeClustersPerDaughter
          ? 1. / static_cast<double>(fillableClusters)
          : 1.;

    const std::size_t size = daughterVectorSize(daughter);

    for (std::size_t index = 0; index < size; ++index)
    {
      const int layer =
          static_cast<int>(daughter.layer->at(index));
      const int radialModule =
          radialModuleFromLayer(layer);
      const int localLayer =
          localLayerBin(layer, radialModule);

      if (radialModule < 0 || localLayer < 0) continue;

      const double x = daughter.clusterX->at(index);
      const double y = daughter.clusterY->at(index);
      const double z = daughter.clusterZ->at(index);
      const double radius = daughter.clusterR->at(index);
      const double phi = wrapPhi(daughter.clusterPhi->at(index));

      if (!std::isfinite(x) ||
          !std::isfinite(y) ||
          !std::isfinite(z) ||
          !std::isfinite(radius) ||
          !std::isfinite(phi) ||
          radius <= 0. ||
          std::abs(z) >= maximumAbsClusterZ)
      {
        continue;
      }

      const PhiBin phiBin =
          determinePhiBin(phi, sectorPhiOffset);
      const int module =
          hardwareModule(phiBin.sector, radialModule);

      if (module < 0 || phiBin.localPhi < 0) continue;

      const auto corrected =
          closestPointOnCircle(correctedCircle, x, y);

      if (!std::isfinite(corrected.first) ||
          !std::isfinite(corrected.second))
      {
        continue;
      }

      const double correctionX = corrected.first - x;
      const double correctionY = corrected.second - y;

      const double radialX = std::cos(phi);
      const double radialY = std::sin(phi);
      const double tangentialX = -std::sin(phi);
      const double tangentialY = std::cos(phi);

      const double deltaR =
          correctionX * radialX +
          correctionY * radialY;
      const double deltaRPhi =
          correctionX * tangentialX +
          correctionY * tangentialY;

      if (!std::isfinite(deltaR) ||
          !std::isfinite(deltaRPhi) ||
          std::abs(deltaR) > maximumCorrection ||
          std::abs(deltaRPhi) > maximumCorrection)
      {
        continue;
      }

      const std::string path = histogramPath(
        selectionName,
        daughter.side,
        daughter.charge,
        iPt,
        module,
        phiBin.localPhi,
        localLayer);

      TH2D*& histogram = voteHistograms[path];
      if (!histogram)
      {
        TDirectory* directory =
            makeDirectoryPath(output, path);

        histogram = bookVoteHistogram(
          directory,
          path,
          voteBins,
          maximumCorrection);
      }

      const double normalDistance =
          std::hypot(correctionX, correctionY);

      if (!(normalDistance > 1e-12) ||
          !std::isfinite(normalDistance))
      {
        continue;
      }

      const double normalX =
          correctionX / normalDistance;
      const double normalY =
          correctionY / normalDistance;

      const double normalR =
          normalX * radialX +
          normalY * radialY;
      const double normalPhi =
          normalX * tangentialX +
          normalY * tangentialY;

      fillAllowedCorrectionLine(
        histogram,
        normalR,
        normalPhi,
        normalDistance,
        maximumCorrection,
        lineSamples,
        weight);

      ++voteCounts[path];
      ++selectedClusters;
    }
  }
}

void MakeK0s2DClusterCorrectionVotes_v3(
    const char* inputDirectory = ".",
    const char* filePattern = "*.root",
    const char* outputName = "k0s_2d_cluster_votes_exact_xy_cut03.root",
    const char* treeName = "pairTree",
    const double beamX = 0.158,
    const double beamY = 0.285,
    const double beamZ = 0.0,
    const double magneticFieldTesla = 1.4,
    const double anchorRadiusMinimum = 46.0,
    const double anchorRadiusMaximum = 54.0,
    const double preferredAnchorRadius = 50.0,
    const char* selectionName = "cut03_baseline",
    const double daughterPtMaximum = 1.50,
    const double kshortPtMaximum = 2.0,
    const double optionalOpeningAngleMinimum = -1.0,
    const double massMinimum = 0.47,
    const double massMaximum = 0.53,
    const int alternatingIterations = 3,
    const double minimumMomentumScale = 0.70,
    const double maximumMomentumScale = 1.30,
    const double maximumAbsRotation = 0.05,
    const double maximumAbsClusterZ = 105.0,
    const double sectorPhiOffset = -TMath::Pi(),
    const double maximumCorrection = 0.8,
    const int voteBins = 121,
    const int lineSamples = 241,
    const bool normalizeClustersPerDaughter = true,
    const Long64_t maximumEntries = -1)
{
  const TString chainPattern =
      TString::Format("%s/%s", inputDirectory, filePattern);

  TChain chain(treeName);
  const int filesAdded = chain.Add(chainPattern);

  if (filesAdded <= 0)
  {
    std::cerr
        << "ERROR: no files matched "
        << chainPattern << std::endl;
    return;
  }

  // Pair-level branches.
  UInt_t candidateMask = 0;
  UChar_t hasKshortDetails = 0;

  Float_t massKshort = 0.F;
  Float_t alpha = 0.F;
  Float_t pairDca = 0.F;
  Float_t pairDedx1 = 0.F;
  Float_t pairDedx2 = 0.F;
  Float_t quality1 = 0.F;
  Float_t quality2 = 0.F;
  Float_t charge1 = 0.F;
  Float_t charge2 = 0.F;

  Float_t pcaX = 0.F;
  Float_t pcaY = 0.F;
  Float_t pcaZ = 0.F;
  Float_t pca1X = 0.F;
  Float_t pca1Y = 0.F;
  Float_t pca1Z = 0.F;
  Float_t pca2X = 0.F;
  Float_t pca2Y = 0.F;
  Float_t pca2Z = 0.F;

  Float_t pairPx1 = 0.F;
  Float_t pairPy1 = 0.F;
  Float_t pairPz1 = 0.F;
  Float_t pairPx2 = 0.F;
  Float_t pairPy2 = 0.F;
  Float_t pairPz2 = 0.F;

  Short_t pairNpoints1 = 0;
  Short_t pairNpoints2 = 0;

  Daughter daughter1;
  Daughter daughter2;

  const std::vector<const char*> requiredBranches = {
    "candidate_mask",
    "has_kshort_daughter_details",
    "mass_Kshort",
    "alpha",
    "pairDCA",
    "dedx_1",
    "dedx_2",
    "quality1",
    "quality2",
    "charge1",
    "charge2",
    "pca_x",
    "pca_y",
    "pca_z",
    "pca1_x",
    "pca1_y",
    "pca1_z",
    "pca2_x",
    "pca2_y",
    "pca2_z",
    "px1",
    "py1",
    "pz1",
    "px2",
    "py2",
    "pz2",
    "npoints1",
    "npoints2"
  };

  for (const char* branchName : requiredBranches)
  {
    if (!chain.GetBranch(branchName))
    {
      std::cerr
          << "ERROR: missing required branch "
          << branchName << std::endl;
      return;
    }
  }

  chain.SetBranchAddress("candidate_mask", &candidateMask);
  chain.SetBranchAddress(
    "has_kshort_daughter_details",
    &hasKshortDetails);

  chain.SetBranchAddress("mass_Kshort", &massKshort);
  chain.SetBranchAddress("alpha", &alpha);
  chain.SetBranchAddress("pairDCA", &pairDca);
  chain.SetBranchAddress("dedx_1", &pairDedx1);
  chain.SetBranchAddress("dedx_2", &pairDedx2);
  chain.SetBranchAddress("quality1", &quality1);
  chain.SetBranchAddress("quality2", &quality2);
  chain.SetBranchAddress("charge1", &charge1);
  chain.SetBranchAddress("charge2", &charge2);

  chain.SetBranchAddress("pca_x", &pcaX);
  chain.SetBranchAddress("pca_y", &pcaY);
  chain.SetBranchAddress("pca_z", &pcaZ);
  chain.SetBranchAddress("pca1_x", &pca1X);
  chain.SetBranchAddress("pca1_y", &pca1Y);
  chain.SetBranchAddress("pca1_z", &pca1Z);
  chain.SetBranchAddress("pca2_x", &pca2X);
  chain.SetBranchAddress("pca2_y", &pca2Y);
  chain.SetBranchAddress("pca2_z", &pca2Z);

  chain.SetBranchAddress("px1", &pairPx1);
  chain.SetBranchAddress("py1", &pairPy1);
  chain.SetBranchAddress("pz1", &pairPz1);
  chain.SetBranchAddress("px2", &pairPx2);
  chain.SetBranchAddress("py2", &pairPy2);
  chain.SetBranchAddress("pz2", &pairPz2);

  chain.SetBranchAddress("npoints1", &pairNpoints1);
  chain.SetBranchAddress("npoints2", &pairNpoints2);

  bindDaughter(chain, "daughter1", daughter1);
  bindDaughter(chain, "daughter2", daughter2);

  std::unique_ptr<TFile> output(
      TFile::Open(outputName, "RECREATE"));

  if (!output || output->IsZombie())
  {
    std::cerr
        << "ERROR: could not create "
        << outputName << std::endl;
    return;
  }

  output->cd();

  const V0Cut selectedCut =
      makeRequestedCut(selectionName ? selectionName : "cut03_baseline");

  std::cout << "Using KShort selection: "
            << selectedCut.name << std::endl;

  // ----------------------------------------------------------
  // Global pair-constraint QA.
  // ----------------------------------------------------------
  TH1D hMassBefore(
    "h_mass_before",
    "K^{0}_{S} mass before constraint;"
    "m_{#pi#pi} [GeV/c^{2}];candidates",
    240, 0.44, 0.56);

  TH1D hMassAfter(
    "h_mass_after",
    "K^{0}_{S} mass after constraint;"
    "m_{#pi#pi} [GeV/c^{2}];candidates",
    240, 0.44, 0.56);

  TH1D hDiraXYBefore(
    "h_dira_xy_before",
    "Transverse DIRA before constraint;DIRA_{xy};candidates",
    220, 0.78, 1.0);

  TH1D hDiraXYAfter(
    "h_dira_xy_after",
    "Transverse DIRA after constraint;DIRA_{xy};candidates",
    220, 0.999, 1.00001);

  TH1D hDira3DBefore(
    "h_dira_3d_before",
    "3D DIRA before constraint;DIRA_{3D};candidates",
    220, 0.78, 1.0);

  TH1D hDira3DAfter(
    "h_dira_3d_after",
    "3D DIRA after constraint;DIRA_{3D};candidates",
    220, 0.78, 1.0);

  TH1D hPairDcaXYBefore(
    "h_pair_dca_xy_before",
    "Transverse circle pair DCA before;"
    "pairDCA_{xy} [cm];candidates",
    240, 0., 2.4);

  TH1D hPairDcaXYAfter(
    "h_pair_dca_xy_after",
    "Transverse circle pair DCA after;"
    "pairDCA_{xy} [cm];candidates",
    240, 0., 2.4);

  TH1D hPairDca3DBefore(
    "h_pair_dca_3d_before",
    "3D line pair DCA before (QA);"
    "pairDCA_{3D} [cm];candidates",
    240, 0., 2.4);

  TH1D hPairDca3DAfter(
    "h_pair_dca_3d_after",
    "3D line pair DCA after (QA);"
    "pairDCA_{3D} [cm];candidates",
    240, 0., 2.4);

  TH1D hScale1(
    "h_momentum_scale1",
    "Daughter 1 momentum scale;s_{1};candidates",
    240, 0.7, 1.3);

  TH1D hScale2(
    "h_momentum_scale2",
    "Daughter 2 momentum scale;s_{2};candidates",
    240, 0.7, 1.3);

  TH2D hScaleCorrelation(
    "h_momentum_scale2_vs_scale1",
    "Momentum-scale correlation;s_{1};s_{2}",
    160, 0.7, 1.3,
    160, 0.7, 1.3);

  TH1D hRotation1(
    "h_rotation1",
    "Daughter 1 transverse rotation;"
    "#delta#phi_{1} [rad];candidates",
    240, -0.06, 0.06);

  TH1D hRotation2(
    "h_rotation2",
    "Daughter 2 transverse rotation;"
    "#delta#phi_{2} [rad];candidates",
    240, -0.06, 0.06);

  TH2D hRotationCorrelation(
    "h_rotation2_vs_rotation1",
    "Rotation correlation;"
    "#delta#phi_{1} [rad];#delta#phi_{2} [rad]",
    160, -0.06, 0.06,
    160, -0.06, 0.06);

  TH1D hOpeningAngle(
    "h_opening_angle",
    "Selected daughter opening angle;"
    "opening angle [rad];candidates",
    180, 0., TMath::Pi());

  TH1D hKshortPt(
    "h_kshort_pt",
    "Selected K^{0}_{S} p_{T};"
    "p_{T}^{K^{0}_{S}} [GeV/c];candidates",
    200, 0., 4.);

  TH1D hDaughterPt(
    "h_daughter_pt",
    "Selected daughter p_{T};"
    "p_{T}^{daughter} [GeV/c];daughters",
    250, 0., 5.);

  TH1D hVoteDeltaR(
    "h_vote_delta_r",
    "All accepted cluster correction votes;"
    "#Delta r [cm];clusters",
    240, -maximumCorrection, maximumCorrection);

  TH1D hVoteDeltaRPhi(
    "h_vote_delta_rphi",
    "All accepted cluster correction votes;"
    "r#Delta#phi [cm];clusters",
    240, -maximumCorrection, maximumCorrection);

  TH1D hVotesPerDaughter(
    "h_votes_per_daughter",
    "Accepted vote clusters per daughter;"
    "clusters;daughters",
    80, 0., 80.);

  TTree fitTree(
    "constraintTree",
    "K0S alternating-constraint fit results");

  Float_t outMassBefore = 0.F;
  Float_t outMassAfter = 0.F;
  Float_t outDiraXYBefore = 0.F;
  Float_t outDiraXYAfter = 0.F;
  Float_t outDira3DBefore = 0.F;
  Float_t outDira3DAfter = 0.F;
  Float_t outPairDcaXYBefore = 0.F;
  Float_t outPairDcaXYAfter = 0.F;
  Float_t outPairDca3DBefore = 0.F;
  Float_t outPairDca3DAfter = 0.F;
  Float_t outScale1 = 1.F;
  Float_t outScale2 = 1.F;
  Float_t outRotation1 = 0.F;
  Float_t outRotation2 = 0.F;
  Float_t outOpeningAngle = 0.F;
  Float_t outPt1 = 0.F;
  Float_t outPt2 = 0.F;
  Float_t outKshortPt = 0.F;
  Int_t outIterations = 0;
  Int_t outValid = 0;

  fitTree.Branch(
    "mass_before",
    &outMassBefore,
    "mass_before/F");
  fitTree.Branch(
    "mass_after",
    &outMassAfter,
    "mass_after/F");
  fitTree.Branch(
    "dira_xy_before",
    &outDiraXYBefore,
    "dira_xy_before/F");
  fitTree.Branch(
    "dira_xy_after",
    &outDiraXYAfter,
    "dira_xy_after/F");
  fitTree.Branch(
    "dira_3d_before",
    &outDira3DBefore,
    "dira_3d_before/F");
  fitTree.Branch(
    "dira_3d_after",
    &outDira3DAfter,
    "dira_3d_after/F");
  fitTree.Branch(
    "pair_dca_xy_before",
    &outPairDcaXYBefore,
    "pair_dca_xy_before/F");
  fitTree.Branch(
    "pair_dca_xy_after",
    &outPairDcaXYAfter,
    "pair_dca_xy_after/F");
  fitTree.Branch(
    "pair_dca_3d_before",
    &outPairDca3DBefore,
    "pair_dca_3d_before/F");
  fitTree.Branch(
    "pair_dca_3d_after",
    &outPairDca3DAfter,
    "pair_dca_3d_after/F");
  fitTree.Branch(
    "momentum_scale1",
    &outScale1,
    "momentum_scale1/F");
  fitTree.Branch(
    "momentum_scale2",
    &outScale2,
    "momentum_scale2/F");
  fitTree.Branch(
    "rotation1",
    &outRotation1,
    "rotation1/F");
  fitTree.Branch(
    "rotation2",
    &outRotation2,
    "rotation2/F");
  fitTree.Branch(
    "opening_angle",
    &outOpeningAngle,
    "opening_angle/F");
  fitTree.Branch("pt1", &outPt1, "pt1/F");
  fitTree.Branch("pt2", &outPt2, "pt2/F");
  fitTree.Branch(
    "kshort_pt",
    &outKshortPt,
    "kshort_pt/F");
  fitTree.Branch(
    "iterations",
    &outIterations,
    "iterations/I");
  fitTree.Branch("valid", &outValid, "valid/I");

  std::map<std::string, TH2D*> voteHistograms;
  std::map<std::string, Long64_t> voteCounts;

  const Long64_t totalEntries = chain.GetEntries();
  const Long64_t entriesToRun =
      maximumEntries > 0
        ? std::min(maximumEntries, totalEntries)
        : totalEntries;

  Long64_t selectedCandidates = 0;
  Long64_t solvedCandidates = 0;
  Long64_t selectedClusters = 0;

  for (Long64_t entry = 0;
       entry < entriesToRun;
       ++entry)
  {
    chain.GetEntry(entry);

    if (entry % 100000 == 0)
    {
      std::cout
          << "Processing "
          << entry << " / "
          << entriesToRun << std::endl;
    }

    if ((candidateMask & 1U) == 0U ||
        hasKshortDetails == 0 ||
        charge1 * charge2 >= 0.F ||
        !daughterVectorsValid(daughter1) ||
        !daughterVectorsValid(daughter2))
    {
      continue;
    }

    const Vec3 momentum1{
      pairPx1,
      pairPy1,
      pairPz1
    };

    const Vec3 momentum2{
      pairPx2,
      pairPy2,
      pairPz2
    };

    const Vec3 secondaryVertex{
      pcaX,
      pcaY,
      pcaZ
    };

    const Vec3 primaryVertex{
      beamX,
      beamY,
      beamZ
    };

    const Vec3 flight =
        subtract(secondaryVertex, primaryVertex);

    const Vec3 daughterPca1{
      pca1X,
      pca1Y,
      pca1Z
    };

    const Vec3 daughterPca2{
      pca2X,
      pca2Y,
      pca2Z
    };

    const double pt1 =
        transverseMomentum(momentum1);
    const double pt2 =
        transverseMomentum(momentum2);
    const Vec3 kshortMomentum =
        add(momentum1, momentum2);
    const double kshortPt =
        transverseMomentum(kshortMomentum);
    const double currentOpeningAngle =
        openingAngle(momentum1, momentum2);
    const double inputDira =
        calculateDira(kshortMomentum, flight);
    const double decayRadius =
        std::hypot(pcaX, pcaY);
    const double deltaPcaZ =
        std::abs(pca1Z - pca2Z);

    const double phiPositive =
        charge1 > 0.F
          ? std::atan2(pairPy1, pairPx1)
          : std::atan2(pairPy2, pairPx2);

    const double phiNegative =
        charge1 > 0.F
          ? std::atan2(pairPy2, pairPx2)
          : std::atan2(pairPy1, pairPx1);

    const double signedDeltaPhi =
        wrapPhi(phiPositive - phiNegative);

    const double minimumSignedDeltaPhi =
        0.8 - 0.4 * std::min(kshortPt, 2.0);

    const bool passV0DeltaPhi =
        signedDeltaPhi >= minimumSignedDeltaPhi;

    const bool passPionDedx =
        std::isfinite(pairDedx1) &&
        std::isfinite(pairDedx2) &&
        pairDedx1 < 400.0 &&
        pairDedx2 < 400.0;

    const bool passOptionalOpeningAngle =
        optionalOpeningAngleMinimum < 0.0 ||
        currentOpeningAngle >= optionalOpeningAngleMinimum;

    const bool passSelection =
        finite(momentum1) &&
        finite(momentum2) &&
        std::isfinite(currentOpeningAngle) &&
        massKshort >= massMinimum &&
        massKshort <= massMaximum &&
        pt1 >= selectedCut.minimumDaughterPt &&
        pt2 >= selectedCut.minimumDaughterPt &&
        pt1 <= daughterPtMaximum &&
        pt2 <= daughterPtMaximum &&
        kshortPt <= kshortPtMaximum &&
        passOptionalOpeningAngle &&
        passV0DeltaPhi &&
        passPionDedx &&
        inputDira >= selectedCut.minimumDira &&
        std::abs(pairDca) <= selectedCut.maximumPairDca &&
        std::abs(pcaZ) < selectedCut.maximumAbsPcaZ &&
        deltaPcaZ < selectedCut.maximumAbsDeltaPcaZ &&
        decayRadius > selectedCut.minimumDecayRadius &&
        std::abs(alpha) < selectedCut.maximumAbsAlpha &&
        quality1 < selectedCut.maximumQuality &&
        quality2 < selectedCut.maximumQuality &&
        pairNpoints1 > selectedCut.minimumNpoints &&
        pairNpoints2 > selectedCut.minimumNpoints;

    if (!passSelection) continue;

    const Anchor anchor1 =
        chooseR2Anchor(
          daughter1,
          anchorRadiusMinimum,
          anchorRadiusMaximum,
          preferredAnchorRadius);

    const Anchor anchor2 =
        chooseR2Anchor(
          daughter2,
          anchorRadiusMinimum,
          anchorRadiusMaximum,
          preferredAnchorRadius);

    if (!anchor1.valid || !anchor2.valid) continue;

    ++selectedCandidates;

    const FitResult fit = solveAlternating(
      momentum1,
      momentum2,
      flight,
      daughterPca1,
      daughterPca2,
      anchor1,
      anchor2,
      daughter1.charge,
      daughter2.charge,
      magneticFieldTesla,
      alternatingIterations,
      minimumMomentumScale,
      maximumMomentumScale,
      maximumAbsRotation);

    outMassBefore = fit.massBefore;
    outMassAfter = fit.massAfter;
    outDiraXYBefore = fit.diraXYBefore;
    outDiraXYAfter = fit.diraXYAfter;
    outDira3DBefore = fit.dira3DBefore;
    outDira3DAfter = fit.dira3DAfter;
    outPairDcaXYBefore = fit.pairDcaXYBefore;
    outPairDcaXYAfter = fit.pairDcaXYAfter;
    outPairDca3DBefore = fit.pairDca3DBefore;
    outPairDca3DAfter = fit.pairDca3DAfter;
    outScale1 = fit.scale1;
    outScale2 = fit.scale2;
    outRotation1 = fit.rotation1;
    outRotation2 = fit.rotation2;
    outOpeningAngle = currentOpeningAngle;
    outPt1 = pt1;
    outPt2 = pt2;
    outKshortPt = kshortPt;
    outIterations = fit.iterations;
    outValid = fit.valid ? 1 : 0;

    fitTree.Fill();

    if (!fit.valid) continue;

    const Circle correctedCircle1 = fit.circle1;
    const Circle correctedCircle2 = fit.circle2;

    if (!correctedCircle1.valid ||
        !correctedCircle2.valid)
    {
      continue;
    }

    hMassBefore.Fill(fit.massBefore);
    hMassAfter.Fill(fit.massAfter);
    hDiraXYBefore.Fill(fit.diraXYBefore);
    hDiraXYAfter.Fill(fit.diraXYAfter);
    hDira3DBefore.Fill(fit.dira3DBefore);
    hDira3DAfter.Fill(fit.dira3DAfter);
    hPairDcaXYBefore.Fill(fit.pairDcaXYBefore);
    hPairDcaXYAfter.Fill(fit.pairDcaXYAfter);
    hPairDca3DBefore.Fill(fit.pairDca3DBefore);
    hPairDca3DAfter.Fill(fit.pairDca3DAfter);
    hScale1.Fill(fit.scale1);
    hScale2.Fill(fit.scale2);
    hScaleCorrelation.Fill(
      fit.scale1,
      fit.scale2);
    hRotation1.Fill(fit.rotation1);
    hRotation2.Fill(fit.rotation2);
    hRotationCorrelation.Fill(
      fit.rotation1,
      fit.rotation2);
    hOpeningAngle.Fill(currentOpeningAngle);
    hKshortPt.Fill(kshortPt);
    hDaughterPt.Fill(pt1);
    hDaughterPt.Fill(pt2);

    const std::size_t votes1 =
        countFillableClusters(
          daughter1,
          maximumAbsClusterZ,
          maximumCorrection,
          correctedCircle1,
          sectorPhiOffset);

    const std::size_t votes2 =
        countFillableClusters(
          daughter2,
          maximumAbsClusterZ,
          maximumCorrection,
          correctedCircle2,
          sectorPhiOffset);

    hVotesPerDaughter.Fill(votes1);
    hVotesPerDaughter.Fill(votes2);

    const Long64_t clustersBefore =
        selectedClusters;

    fillDaughterVotes(
      selectedCut.name,
      daughter1,
      pt1,
      correctedCircle1,
      output.get(),
      voteHistograms,
      voteCounts,
      maximumAbsClusterZ,
      maximumCorrection,
      voteBins,
      lineSamples,
      sectorPhiOffset,
      normalizeClustersPerDaughter,
      selectedClusters);

    fillDaughterVotes(
      selectedCut.name,
      daughter2,
      pt2,
      correctedCircle2,
      output.get(),
      voteHistograms,
      voteCounts,
      maximumAbsClusterZ,
      maximumCorrection,
      voteBins,
      lineSamples,
      sectorPhiOffset,
      normalizeClustersPerDaughter,
      selectedClusters);

    if (selectedClusters > clustersBefore)
      ++solvedCandidates;
  }

  // Fill global vote projections from every final histogram.
  for (const auto& [path, histogram] : voteHistograms)
  {
    if (!histogram) continue;

    for (int binX = 1;
         binX <= histogram->GetNbinsX();
         ++binX)
    {
      for (int binY = 1;
           binY <= histogram->GetNbinsY();
           ++binY)
      {
        const double weight =
            histogram->GetBinContent(binX, binY);

        if (!(weight > 0.)) continue;

        hVoteDeltaR.Fill(
          histogram->GetXaxis()->GetBinCenter(binX),
          weight);

        hVoteDeltaRPhi.Fill(
          histogram->GetYaxis()->GetBinCenter(binY),
          weight);
      }
    }
  }

  output->cd();

  TNamed constraintMethod(
    "constraint_method",
    "Exact transverse pointing plus exact PDG mass for the two momentum "
    "scales; transverse circle pairDCA for the two rotations; full allowed "
    "cluster correction lines are filled.");
  constraintMethod.Write();

  TNamed selectionConvention(
    "kshort_selection",
    TString::Format(
      "%s; dedx_1<400; dedx_2<400; "
      "signed DeltaPhi >= 0.8 - 0.4*min(KShort pT,2); "
      "optional opening-angle threshold = %.6g rad",
      selectedCut.name.c_str(),
      optionalOpeningAngleMinimum).Data());
  selectionConvention.Write();

  TNamed moduleConvention(
    "module_numbering",
    "modules 0-11 are R1 sectors 0-11; "
    "modules 12-23 are R2 sectors 0-11; "
    "modules 24-35 are R3 sectors 0-11");
  moduleConvention.Write();

  TNamed weightingConvention(
    "weighting_convention",
    normalizeClustersPerDaughter
      ? "Each daughter contributes unit total weight distributed over its "
        "accepted clusters. Final notebook additionally normalizes every "
        "charge-pT histogram before combining categories."
      : "Every accepted cluster contributes unit weight. Final notebook "
        "normalizes every charge-pT histogram before combining categories.");
  weightingConvention.Write();

  output->Write();
  output->Close();

  std::cout
      << "Input entries: " << totalEntries << '\n'
      << "Selected candidates: " << selectedCandidates << '\n'
      << "Solved candidates with votes: " << solvedCandidates << '\n'
      << "Accepted cluster votes: " << selectedClusters << '\n'
      << "Booked vote histograms: " << voteHistograms.size() << '\n'
      << "Wrote: " << outputName << std::endl;
}
