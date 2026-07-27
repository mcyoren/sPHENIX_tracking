// MakeK0s2DClusterCorrectionVotes.C
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
// root -l -b -q 'MakeK0s2DClusterCorrectionVotes.C(
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
    double diraBefore = -2.;
    double diraAfter = -2.;
    double pairDcaBefore = 0.;
    double pairDcaAfter = 0.;
    Vec3 momentum1;
    Vec3 momentum2;
  };

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
      const int side,
      const int charge,
      const int iPt,
      const int module,
      const int localPhi,
      const int localLayer)
  {
    return TString::Format(
      "side%d/%s/%s/module_%02d/phi_%d_layer_%d",
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

  double objective(
      const Vec3& momentum1,
      const Vec3& momentum2,
      const Vec3& flight,
      const Vec3& pca1,
      const Vec3& pca2,
      const double massSigma,
      const double pointingSigma,
      const double pairDcaSigma,
      const double scale1,
      const double scale2,
      const double rotation1,
      const double rotation2,
      const double scaleRegularization,
      const double rotationRegularization)
  {
    const double mass =
        invariantMass(momentum1, momentum2);
    const double currentDira =
        calculateDira(add(momentum1, momentum2), flight);
    const double currentPairDca =
        lineLineDca(pca1, momentum1, pca2, momentum2);

    if (!std::isfinite(mass) ||
        !std::isfinite(currentDira) ||
        !std::isfinite(currentPairDca))
    {
      return std::numeric_limits<double>::max();
    }

    const double massTerm =
        (mass - kKshortMass) / massSigma;
    const double pointingTerm =
        (1. - currentDira) / pointingSigma;
    const double pairDcaTerm =
        currentPairDca / pairDcaSigma;

    return massTerm * massTerm +
           pointingTerm * pointingTerm +
           pairDcaTerm * pairDcaTerm +
           scaleRegularization *
             ((scale1 - 1.) * (scale1 - 1.) +
              (scale2 - 1.) * (scale2 - 1.)) +
           rotationRegularization *
             (rotation1 * rotation1 +
              rotation2 * rotation2);
  }

  FitResult solveAlternating(
      const Vec3& originalMomentum1,
      const Vec3& originalMomentum2,
      const Vec3& flight,
      const Vec3& pca1,
      const Vec3& pca2,
      const int requestedIterations,
      const double minimumScale,
      const double maximumScale,
      const double maximumRotation,
      const double massSigma,
      const double pointingSigma,
      const double pairDcaSigma)
  {
    FitResult result;

    result.massBefore =
        invariantMass(originalMomentum1, originalMomentum2);
    result.diraBefore =
        calculateDira(add(originalMomentum1, originalMomentum2), flight);
    result.pairDcaBefore =
        lineLineDca(
          pca1, originalMomentum1,
          pca2, originalMomentum2);

    double scale1 = 1.;
    double scale2 = 1.;
    double rotation1 = 0.;
    double rotation2 = 0.;

    // A modest regularization prevents weakly constrained fits from running
    // to scale or rotation boundaries.
    constexpr double scaleRegularization = 0.10;
    constexpr double rotationRegularization = 25.0;

    for (int iteration = 0;
         iteration < requestedIterations;
         ++iteration)
    {
      // ------------------------------------------------------
      // Step A: fixed rotations, solve the two momentum scales.
      // ------------------------------------------------------
      double centerScale1 = scale1;
      double centerScale2 = scale2;
      double halfScale =
          iteration == 0
            ? 0.20
            : 0.05;

      for (int refinement = 0; refinement < 4; ++refinement)
      {
        double bestValue =
            std::numeric_limits<double>::max();
        double bestScale1 = centerScale1;
        double bestScale2 = centerScale2;

        constexpr int steps = 35;

        for (int index1 = 0; index1 <= steps; ++index1)
        {
          const double fraction1 =
              static_cast<double>(index1) / steps;
          const double candidateScale1 =
              std::clamp(
                centerScale1 - halfScale +
                  2. * halfScale * fraction1,
                minimumScale,
                maximumScale);

          for (int index2 = 0; index2 <= steps; ++index2)
          {
            const double fraction2 =
                static_cast<double>(index2) / steps;
            const double candidateScale2 =
                std::clamp(
                  centerScale2 - halfScale +
                    2. * halfScale * fraction2,
                  minimumScale,
                  maximumScale);

            const Vec3 candidateMomentum1 =
                rotateZ(
                  multiply(originalMomentum1, candidateScale1),
                  rotation1);
            const Vec3 candidateMomentum2 =
                rotateZ(
                  multiply(originalMomentum2, candidateScale2),
                  rotation2);

            const double value = objective(
              candidateMomentum1,
              candidateMomentum2,
              flight,
              pca1,
              pca2,
              massSigma,
              pointingSigma,
              pairDcaSigma,
              candidateScale1,
              candidateScale2,
              rotation1,
              rotation2,
              scaleRegularization,
              rotationRegularization);

            if (value < bestValue)
            {
              bestValue = value;
              bestScale1 = candidateScale1;
              bestScale2 = candidateScale2;
            }
          }
        }

        centerScale1 = bestScale1;
        centerScale2 = bestScale2;
        halfScale *= 0.20;
      }

      scale1 = centerScale1;
      scale2 = centerScale2;

      // ------------------------------------------------------
      // Step B: fixed momentum scales, solve small rotations.
      // ------------------------------------------------------
      double centerRotation1 = rotation1;
      double centerRotation2 = rotation2;
      double halfRotation =
          iteration == 0
            ? maximumRotation
            : 0.25 * maximumRotation;

      for (int refinement = 0; refinement < 4; ++refinement)
      {
        double bestValue =
            std::numeric_limits<double>::max();
        double bestRotation1 = centerRotation1;
        double bestRotation2 = centerRotation2;

        constexpr int steps = 35;

        for (int index1 = 0; index1 <= steps; ++index1)
        {
          const double fraction1 =
              static_cast<double>(index1) / steps;
          const double candidateRotation1 =
              std::clamp(
                centerRotation1 - halfRotation +
                  2. * halfRotation * fraction1,
                -maximumRotation,
                maximumRotation);

          for (int index2 = 0; index2 <= steps; ++index2)
          {
            const double fraction2 =
                static_cast<double>(index2) / steps;
            const double candidateRotation2 =
                std::clamp(
                  centerRotation2 - halfRotation +
                    2. * halfRotation * fraction2,
                  -maximumRotation,
                  maximumRotation);

            const Vec3 candidateMomentum1 =
                rotateZ(
                  multiply(originalMomentum1, scale1),
                  candidateRotation1);
            const Vec3 candidateMomentum2 =
                rotateZ(
                  multiply(originalMomentum2, scale2),
                  candidateRotation2);

            const double value = objective(
              candidateMomentum1,
              candidateMomentum2,
              flight,
              pca1,
              pca2,
              massSigma,
              pointingSigma,
              pairDcaSigma,
              scale1,
              scale2,
              candidateRotation1,
              candidateRotation2,
              scaleRegularization,
              rotationRegularization);

            if (value < bestValue)
            {
              bestValue = value;
              bestRotation1 = candidateRotation1;
              bestRotation2 = candidateRotation2;
            }
          }
        }

        centerRotation1 = bestRotation1;
        centerRotation2 = bestRotation2;
        halfRotation *= 0.20;
      }

      rotation1 = centerRotation1;
      rotation2 = centerRotation2;
      result.iterations = iteration + 1;
    }

    result.scale1 = scale1;
    result.scale2 = scale2;
    result.rotation1 = rotation1;
    result.rotation2 = rotation2;

    result.momentum1 =
        rotateZ(
          multiply(originalMomentum1, scale1),
          rotation1);
    result.momentum2 =
        rotateZ(
          multiply(originalMomentum2, scale2),
          rotation2);

    result.massAfter =
        invariantMass(result.momentum1, result.momentum2);
    result.diraAfter =
        calculateDira(add(result.momentum1, result.momentum2), flight);
    result.pairDcaAfter =
        lineLineDca(
          pca1, result.momentum1,
          pca2, result.momentum2);

    result.valid =
        finite(result.momentum1) &&
        finite(result.momentum2) &&
        std::isfinite(result.massAfter) &&
        std::isfinite(result.diraAfter) &&
        std::isfinite(result.pairDcaAfter) &&
        scale1 >= minimumScale &&
        scale1 <= maximumScale &&
        scale2 >= minimumScale &&
        scale2 <= maximumScale &&
        std::abs(rotation1) <= maximumRotation &&
        std::abs(rotation2) <= maximumRotation;

    return result;
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
      const Daughter& daughter,
      const double originalPt,
      const Circle& correctedCircle,
      TFile* output,
      std::map<std::string, TH2D*>& voteHistograms,
      std::map<std::string, Long64_t>& voteCounts,
      const double maximumAbsClusterZ,
      const double maximumCorrection,
      const int voteBins,
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

      histogram->Fill(deltaR, deltaRPhi, weight);
      ++voteCounts[path];
      ++selectedClusters;
    }
  }
}

void MakeK0s2DClusterCorrectionVotes(
    const char* inputDirectory = ".",
    const char* filePattern = "*.root",
    const char* outputName = "k0s_2d_cluster_votes.root",
    const char* treeName = "pairTree",
    const double beamX = 0.158,
    const double beamY = 0.285,
    const double beamZ = 0.0,
    const double magneticFieldTesla = 1.4,
    const double anchorRadiusMinimum = 46.0,
    const double anchorRadiusMaximum = 54.0,
    const double preferredAnchorRadius = 50.0,
    const double daughterPtMinimum = 0.20,
    const double daughterPtMaximum = 1.50,
    const double kshortPtMaximum = 2.0,
    const double openingAngleMinimum = 0.50,
    const double massMinimum = 0.47,
    const double massMaximum = 0.53,
    const double minimumInputDira = 0.85,
    const double maximumInputPairDca = 2.0,
    const double maximumAbsPcaZ = 15.0,
    const double maximumAbsDeltaPcaZ = 0.50,
    const double minimumDecayRadius = 2.0,
    const double maximumAbsAlpha = 0.99,
    const double maximumQuality = 15.0,
    const int minimumNpoints = 30,
    const int alternatingIterations = 3,
    const double minimumMomentumScale = 0.70,
    const double maximumMomentumScale = 1.30,
    const double maximumAbsRotation = 0.05,
    const double massConstraintSigma = 0.002,
    const double pointingConstraintSigma = 0.002,
    const double pairDcaConstraintSigma = 0.05,
    const double maximumAbsClusterZ = 105.0,
    const double sectorPhiOffset = -TMath::Pi(),
    const double maximumCorrection = 0.8,
    const int voteBins = 121,
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

  TH1D hDiraBefore(
    "h_dira_before",
    "DIRA before constraint;DIRA;candidates",
    220, 0.78, 1.0);

  TH1D hDiraAfter(
    "h_dira_after",
    "DIRA after constraint;DIRA;candidates",
    220, 0.78, 1.0);

  TH1D hPairDcaBefore(
    "h_pair_dca_before",
    "Pair DCA before constraint;"
    "pair DCA [cm];candidates",
    240, 0., 2.4);

  TH1D hPairDcaAfter(
    "h_pair_dca_after",
    "Pair DCA after constraint;"
    "pair DCA [cm];candidates",
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
  Float_t outDiraBefore = 0.F;
  Float_t outDiraAfter = 0.F;
  Float_t outPairDcaBefore = 0.F;
  Float_t outPairDcaAfter = 0.F;
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
    "dira_before",
    &outDiraBefore,
    "dira_before/F");
  fitTree.Branch(
    "dira_after",
    &outDiraAfter,
    "dira_after/F");
  fitTree.Branch(
    "pair_dca_before",
    &outPairDcaBefore,
    "pair_dca_before/F");
  fitTree.Branch(
    "pair_dca_after",
    &outPairDcaAfter,
    "pair_dca_after/F");
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

    const bool passSelection =
        finite(momentum1) &&
        finite(momentum2) &&
        std::isfinite(currentOpeningAngle) &&
        massKshort >= massMinimum &&
        massKshort <= massMaximum &&
        pt1 >= daughterPtMinimum &&
        pt2 >= daughterPtMinimum &&
        pt1 <= daughterPtMaximum &&
        pt2 <= daughterPtMaximum &&
        kshortPt <= kshortPtMaximum &&
        currentOpeningAngle >= openingAngleMinimum &&
        inputDira >= minimumInputDira &&
        std::abs(pairDca) <= maximumInputPairDca &&
        std::abs(pcaZ) < maximumAbsPcaZ &&
        deltaPcaZ < maximumAbsDeltaPcaZ &&
        decayRadius > minimumDecayRadius &&
        std::abs(alpha) < maximumAbsAlpha &&
        quality1 < maximumQuality &&
        quality2 < maximumQuality &&
        pairNpoints1 > minimumNpoints &&
        pairNpoints2 > minimumNpoints;

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
      alternatingIterations,
      minimumMomentumScale,
      maximumMomentumScale,
      maximumAbsRotation,
      massConstraintSigma,
      pointingConstraintSigma,
      pairDcaConstraintSigma);

    outMassBefore = fit.massBefore;
    outMassAfter = fit.massAfter;
    outDiraBefore = fit.diraBefore;
    outDiraAfter = fit.diraAfter;
    outPairDcaBefore = fit.pairDcaBefore;
    outPairDcaAfter = fit.pairDcaAfter;
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

    const Circle correctedCircle1 =
        circleFromAnchor(
          anchor1,
          fit.scale1,
          fit.rotation1,
          daughter1.charge,
          magneticFieldTesla);

    const Circle correctedCircle2 =
        circleFromAnchor(
          anchor2,
          fit.scale2,
          fit.rotation2,
          daughter2.charge,
          magneticFieldTesla);

    if (!correctedCircle1.valid ||
        !correctedCircle2.valid)
    {
      continue;
    }

    hMassBefore.Fill(fit.massBefore);
    hMassAfter.Fill(fit.massAfter);
    hDiraBefore.Fill(fit.diraBefore);
    hDiraAfter.Fill(fit.diraAfter);
    hPairDcaBefore.Fill(fit.pairDcaBefore);
    hPairDcaAfter.Fill(fit.pairDcaAfter);
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
      daughter1,
      pt1,
      correctedCircle1,
      output.get(),
      voteHistograms,
      voteCounts,
      maximumAbsClusterZ,
      maximumCorrection,
      voteBins,
      sectorPhiOffset,
      normalizeClustersPerDaughter,
      selectedClusters);

    fillDaughterVotes(
      daughter2,
      pt2,
      correctedCircle2,
      output.get(),
      voteHistograms,
      voteCounts,
      maximumAbsClusterZ,
      maximumCorrection,
      voteBins,
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
