// ExtractK0sMomentumScaleSamples.C
//
// Extract per-daughter momentum-scale measurements from clean
// K0S -> pi+ pi- candidates.
//
// With the daughter directions fixed in phi:
//   1. transverse pointing fixes scale2 / scale1 exactly;
//   2. the PDG K0S mass fixes the remaining common scale.
//
// The output tree contains one row per daughter:
//   side, charge, pt, phi, eta, momentum_scale
//
// These samples are intended for a smooth phase-space map
//
//   scale = scale(side, charge, pt, phi, eta)
//
// built by the companion notebook.
//
// Example:
// root -l -b -q 'ExtractK0sMomentumScaleSamples.C(
//   "/path/to/files","p_v0_*.root",
//   "output/k0s_momentum_scale_samples_cut03.root",
//   "pairTree","cut03_baseline")'

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TTree.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
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

  struct V0Cut
  {
    std::string name;
    double maximumAbsPcaZ = 0.;
    double maximumAbsDeltaPcaZ = 0.;
    double minimumDaughterPt = 0.;
    double minimumDecayRadius = 0.;
    double maximumAbsAlpha = 0.;
    double maximumPairDca = 0.;
    double minimumDira3D = 0.;
    double maximumQuality = 0.;
    int minimumNpoints = 0;
  };

  struct ScaleSolution
  {
    bool valid = false;
    double scale1 = 1.;
    double scale2 = 1.;
    double massAfter = 0.;
    double diraXYAfter = -2.;
    double scaleRatio = 1.;
  };

  V0Cut makeCut(const std::string& name)
  {
    if (name == "cut03_baseline")
    {
      return {
        "cut03_baseline",
        15.0,  // |pair vertex z|
        0.50,  // |daughter PCA z difference|
        0.20,  // daughter pT
        2.0,   // decay radius
        0.99,  // |alpha|
        2.00,  // pairDCA
        0.85,  // 3D DIRA used by established QA selection
        15.0,  // quality
        30     // npoints
      };
    }

    if (name == "cut07_pairDCA_5mm")
    {
      return {
        "cut07_pairDCA_5mm",
        10.0,
        0.20,
        0.20,
        2.0,
        0.99,
        0.50,
        0.95,
        10.0,
        32
      };
    }

    std::cerr
        << "WARNING: unknown cut '" << name
        << "'. Using cut03_baseline." << std::endl;

    return makeCut("cut03_baseline");
  }

  Vec3 add(const Vec3& first, const Vec3& second)
  {
    return {
      first.x + second.x,
      first.y + second.y,
      first.z + second.z
    };
  }

  Vec3 multiply(const Vec3& vector, const double value)
  {
    return {
      value * vector.x,
      value * vector.y,
      value * vector.z
    };
  }

  double norm2(const Vec3& vector)
  {
    return
      vector.x * vector.x +
      vector.y * vector.y +
      vector.z * vector.z;
  }

  double norm(const Vec3& vector)
  {
    return std::sqrt(norm2(vector));
  }

  double transverseMomentum(const Vec3& vector)
  {
    return std::hypot(vector.x, vector.y);
  }

  bool finite(const Vec3& vector)
  {
    return
      std::isfinite(vector.x) &&
      std::isfinite(vector.y) &&
      std::isfinite(vector.z);
  }

  double wrapPhi(double phi)
  {
    while (phi >= TMath::Pi()) phi -= kTwoPi;
    while (phi < -TMath::Pi()) phi += kTwoPi;
    return phi;
  }

  double pseudorapidity(const Vec3& momentum)
  {
    const double magnitude = norm(momentum);
    const double denominator = magnitude - momentum.z;
    const double numerator = magnitude + momentum.z;

    if (!(denominator > 0.) || !(numerator > 0.))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    return 0.5 * std::log(numerator / denominator);
  }

  double invariantMass(
      const Vec3& momentum1,
      const Vec3& momentum2)
  {
    const double energy1 =
        std::sqrt(norm2(momentum1) + kPionMass * kPionMass);
    const double energy2 =
        std::sqrt(norm2(momentum2) + kPionMass * kPionMass);

    const Vec3 total = add(momentum1, momentum2);
    const double mass2 =
        (energy1 + energy2) * (energy1 + energy2) -
        norm2(total);

    return mass2 > 0. ? std::sqrt(mass2) : 0.;
  }

  double dira3D(
      const Vec3& momentum,
      const Vec3& flight)
  {
    const double denominator = norm(momentum) * norm(flight);
    if (!(denominator > 0.)) return -2.;

    const double value =
        (momentum.x * flight.x +
         momentum.y * flight.y +
         momentum.z * flight.z) /
        denominator;

    return std::clamp(value, -1., 1.);
  }

  double diraXY(
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

    const double value =
        (momentum.x * flight.x +
         momentum.y * flight.y) /
        (momentumMagnitude * flightMagnitude);

    return std::clamp(value, -1., 1.);
  }

  double transverseCross(
      const Vec3& first,
      const Vec3& second)
  {
    return first.x * second.y -
           first.y * second.x;
  }

  ScaleSolution solveScales(
      const Vec3& momentum1,
      const Vec3& momentum2,
      const Vec3& flight,
      const double minimumScale,
      const double maximumScale,
      const double massTolerance,
      const double diraXYTolerance)
  {
    ScaleSolution result;

    const double denominator =
        transverseCross(momentum2, flight);

    if (!finite(momentum1) ||
        !finite(momentum2) ||
        std::abs(denominator) < 1e-12)
    {
      return result;
    }

    // Exact transverse pointing:
    //
    //   s1 cross(p1,L) + s2 cross(p2,L) = 0.
    //
    const double ratio =
        -transverseCross(momentum1, flight) /
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

    auto massResidual = [&](const double commonScale)
    {
      const Vec3 corrected1 =
          multiply(momentum1, commonScale);
      const Vec3 corrected2 =
          multiply(momentum2, ratio * commonScale);

      return invariantMass(corrected1, corrected2) -
             kKshortMass;
    };

    constexpr int scanSteps = 600;

    double previousScale = lower;
    double previousValue = massResidual(previousScale);

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
          massResidual(currentScale);

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

    if (!bracketFound)
    {
      return result;
    }

    double low = bracketLow;
    double high = bracketHigh;
    double lowValue = massResidual(low);
    double commonScale = bestScale;

    for (int iteration = 0; iteration < 100; ++iteration)
    {
      const double middle = 0.5 * (low + high);
      const double middleValue = massResidual(middle);

      if (!std::isfinite(middleValue))
      {
        return result;
      }

      commonScale = middle;

      if (std::abs(middleValue) < massTolerance)
      {
        break;
      }

      if (lowValue * middleValue <= 0.)
      {
        high = middle;
      }
      else
      {
        low = middle;
        lowValue = middleValue;
      }
    }

    const Vec3 corrected1 =
        multiply(momentum1, commonScale);
    const Vec3 corrected2 =
        multiply(momentum2, ratio * commonScale);
    const Vec3 total = add(corrected1, corrected2);

    const double forward =
        total.x * flight.x +
        total.y * flight.y;

    result.scale1 = commonScale;
    result.scale2 = ratio * commonScale;
    result.scaleRatio = ratio;
    result.massAfter =
        invariantMass(corrected1, corrected2);
    result.diraXYAfter =
        diraXY(total, flight);

    result.valid =
        std::isfinite(result.scale1) &&
        std::isfinite(result.scale2) &&
        result.scale1 >= minimumScale &&
        result.scale1 <= maximumScale &&
        result.scale2 >= minimumScale &&
        result.scale2 <= maximumScale &&
        forward > 0. &&
        std::abs(result.massAfter - kKshortMass) <
          10. * massTolerance &&
        result.diraXYAfter >
          1. - diraXYTolerance;

    return result;
  }

  int inferSide(
      const int detailedSide,
      const Vec3& momentum)
  {
    if (detailedSide == 0 || detailedSide == 1)
    {
      return detailedSide;
    }

    // Fallback used only if the detailed side is unavailable.
    return momentum.z >= 0. ? 1 : 0;
  }
}

void ExtractK0sMomentumScaleSamples(
    const char* inputDirectory = ".",
    const char* filePattern = "*.root",
    const char* outputName =
      "k0s_momentum_scale_samples_cut03.root",
    const char* treeName = "pairTree",
    const char* selectionName = "cut03_baseline",
    const double beamX = 0.158,
    const double beamY = 0.285,
    const double beamZ = 0.0,
    const double massMinimum = 0.47,
    const double massMaximum = 0.53,
    const double daughterPtMaximum = 5.0,
    const double kshortPtMaximum = 3.0,
    const double optionalOpeningAngleMinimum = -1.0,
    const double minimumMomentumScale = 0.70,
    const double maximumMomentumScale = 1.30,
    const double massTolerance = 1e-9,
    const double diraXYTolerance = 1e-8,
    const Long64_t maximumEntries = -1)
{
  const TString chainPattern =
      TString::Format(
        "%s/%s",
        inputDirectory,
        filePattern);

  TChain chain(treeName);
  const int numberOfFiles = chain.Add(chainPattern);

  if (numberOfFiles <= 0)
  {
    std::cerr
        << "ERROR: no files matched "
        << chainPattern << std::endl;
    return;
  }

  const V0Cut cut =
      makeCut(selectionName ? selectionName : "cut03_baseline");

  UInt_t candidateMask = 0;
  UChar_t hasKshortDetails = 0;

  Float_t massKshort = 0.F;
  Float_t alpha = 0.F;
  Float_t pairDca = 0.F;
  Float_t dedx1 = 0.F;
  Float_t dedx2 = 0.F;
  Float_t quality1 = 0.F;
  Float_t quality2 = 0.F;
  Float_t charge1 = 0.F;
  Float_t charge2 = 0.F;

  Float_t pcaX = 0.F;
  Float_t pcaY = 0.F;
  Float_t pcaZ = 0.F;
  Float_t pca1Z = 0.F;
  Float_t pca2Z = 0.F;

  Float_t px1 = 0.F;
  Float_t py1 = 0.F;
  Float_t pz1 = 0.F;
  Float_t px2 = 0.F;
  Float_t py2 = 0.F;
  Float_t pz2 = 0.F;

  Short_t npoints1 = 0;
  Short_t npoints2 = 0;

  Int_t daughter1Side = -1;
  Int_t daughter2Side = -1;

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
    "pca1_z",
    "pca2_z",
    "px1",
    "py1",
    "pz1",
    "px2",
    "py2",
    "pz2",
    "npoints1",
    "npoints2",
    "daughter1_side",
    "daughter2_side"
  };

  for (const char* branchName : requiredBranches)
  {
    if (!chain.GetBranch(branchName))
    {
      std::cerr
          << "ERROR: missing branch "
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
  chain.SetBranchAddress("dedx_1", &dedx1);
  chain.SetBranchAddress("dedx_2", &dedx2);
  chain.SetBranchAddress("quality1", &quality1);
  chain.SetBranchAddress("quality2", &quality2);
  chain.SetBranchAddress("charge1", &charge1);
  chain.SetBranchAddress("charge2", &charge2);

  chain.SetBranchAddress("pca_x", &pcaX);
  chain.SetBranchAddress("pca_y", &pcaY);
  chain.SetBranchAddress("pca_z", &pcaZ);
  chain.SetBranchAddress("pca1_z", &pca1Z);
  chain.SetBranchAddress("pca2_z", &pca2Z);

  chain.SetBranchAddress("px1", &px1);
  chain.SetBranchAddress("py1", &py1);
  chain.SetBranchAddress("pz1", &pz1);
  chain.SetBranchAddress("px2", &px2);
  chain.SetBranchAddress("py2", &py2);
  chain.SetBranchAddress("pz2", &pz2);

  chain.SetBranchAddress("npoints1", &npoints1);
  chain.SetBranchAddress("npoints2", &npoints2);

  chain.SetBranchAddress(
    "daughter1_side",
    &daughter1Side);
  chain.SetBranchAddress(
    "daughter2_side",
    &daughter2Side);

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

  TH1D hMassBefore(
    "h_mass_before",
    "Selected K^{0}_{S} mass before;"
    "m_{#pi#pi} [GeV/c^{2}];candidates",
    240, 0.44, 0.56);

  TH1D hMassAfter(
    "h_mass_after",
    "K^{0}_{S} mass after exact scale solution;"
    "m_{#pi#pi} [GeV/c^{2}];candidates",
    240, 0.4970, 0.4982);

  TH1D hDiraXYBefore(
    "h_dira_xy_before",
    "Transverse DIRA before;DIRA_{xy};candidates",
    220, 0.78, 1.0);

  TH1D hDiraXYAfter(
    "h_dira_xy_after",
    "Transverse DIRA after exact scale solution;"
    "DIRA_{xy};candidates",
    220, 0.999, 1.00001);

  TH1D hDira3DBefore(
    "h_dira_3d_before",
    "3D DIRA before;DIRA_{3D};candidates",
    220, 0.78, 1.0);

  TH1D hScale1(
    "h_momentum_scale1",
    "Daughter 1 momentum scale;s_{1};candidates",
    240, minimumMomentumScale, maximumMomentumScale);

  TH1D hScale2(
    "h_momentum_scale2",
    "Daughter 2 momentum scale;s_{2};candidates",
    240, minimumMomentumScale, maximumMomentumScale);

  TH2D hScale2VsScale1(
    "h_momentum_scale2_vs_scale1",
    "Momentum-scale correlation;s_{1};s_{2}",
    160, minimumMomentumScale, maximumMomentumScale,
    160, minimumMomentumScale, maximumMomentumScale);

  TH1D hFailureCode(
    "h_failure_code",
    "Scale extraction result;code;candidates",
    5, -0.5, 4.5);

  hFailureCode.GetXaxis()->SetBinLabel(1, "selected");
  hFailureCode.GetXaxis()->SetBinLabel(2, "invalid scale");
  hFailureCode.GetXaxis()->SetBinLabel(3, "non-finite");
  hFailureCode.GetXaxis()->SetBinLabel(4, "side invalid");
  hFailureCode.GetXaxis()->SetBinLabel(5, "written");

  TTree scaleTree(
    "momentumScaleTree",
    "Per-daughter KShort momentum-scale samples");

  Int_t outSide = -1;
  Int_t outCharge = 0;
  Int_t outDaughterIndex = 0;

  Float_t outPt = 0.F;
  Float_t outPhi = 0.F;
  Float_t outEta = 0.F;
  Float_t outP = 0.F;
  Float_t outMomentumScale = 1.F;

  Float_t outOtherPt = 0.F;
  Float_t outOtherMomentumScale = 1.F;

  Float_t outKshortPt = 0.F;
  Float_t outMassBefore = 0.F;
  Float_t outMassAfter = 0.F;
  Float_t outDiraXYBefore = 0.F;
  Float_t outDiraXYAfter = 0.F;
  Float_t outDira3DBefore = 0.F;
  Float_t outPairDca = 0.F;
  Float_t outDecayRadius = 0.F;
  Float_t outAlpha = 0.F;

  scaleTree.Branch("side", &outSide, "side/I");
  scaleTree.Branch("charge", &outCharge, "charge/I");
  scaleTree.Branch(
    "daughter_index",
    &outDaughterIndex,
    "daughter_index/I");

  scaleTree.Branch("pt", &outPt, "pt/F");
  scaleTree.Branch("phi", &outPhi, "phi/F");
  scaleTree.Branch("eta", &outEta, "eta/F");
  scaleTree.Branch("p", &outP, "p/F");
  scaleTree.Branch(
    "momentum_scale",
    &outMomentumScale,
    "momentum_scale/F");

  scaleTree.Branch(
    "other_pt",
    &outOtherPt,
    "other_pt/F");
  scaleTree.Branch(
    "other_momentum_scale",
    &outOtherMomentumScale,
    "other_momentum_scale/F");

  scaleTree.Branch(
    "kshort_pt",
    &outKshortPt,
    "kshort_pt/F");
  scaleTree.Branch(
    "mass_before",
    &outMassBefore,
    "mass_before/F");
  scaleTree.Branch(
    "mass_after",
    &outMassAfter,
    "mass_after/F");
  scaleTree.Branch(
    "dira_xy_before",
    &outDiraXYBefore,
    "dira_xy_before/F");
  scaleTree.Branch(
    "dira_xy_after",
    &outDiraXYAfter,
    "dira_xy_after/F");
  scaleTree.Branch(
    "dira_3d_before",
    &outDira3DBefore,
    "dira_3d_before/F");
  scaleTree.Branch(
    "pair_dca",
    &outPairDca,
    "pair_dca/F");
  scaleTree.Branch(
    "decay_radius",
    &outDecayRadius,
    "decay_radius/F");
  scaleTree.Branch("alpha", &outAlpha, "alpha/F");

  const Long64_t totalEntries = chain.GetEntries();
  const Long64_t entriesToRun =
      maximumEntries > 0
        ? std::min(maximumEntries, totalEntries)
        : totalEntries;

  Long64_t selectedCandidates = 0;
  Long64_t solvedCandidates = 0;
  Long64_t writtenDaughters = 0;

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
        charge1 * charge2 >= 0.F)
    {
      continue;
    }

    const Vec3 momentum1{px1, py1, pz1};
    const Vec3 momentum2{px2, py2, pz2};
    const Vec3 flight{
      pcaX - beamX,
      pcaY - beamY,
      pcaZ - beamZ
    };

    const double pt1 = transverseMomentum(momentum1);
    const double pt2 = transverseMomentum(momentum2);
    const Vec3 totalMomentum = add(momentum1, momentum2);
    const double kshortPt =
        transverseMomentum(totalMomentum);
    const double decayRadius =
        std::hypot(pcaX, pcaY);
    const double opening =
        std::acos(std::clamp(
          (momentum1.x * momentum2.x +
           momentum1.y * momentum2.y +
           momentum1.z * momentum2.z) /
          (norm(momentum1) * norm(momentum2)),
          -1.,
          1.));

    const double phiPositive =
        charge1 > 0.F
          ? std::atan2(py1, px1)
          : std::atan2(py2, px2);
    const double phiNegative =
        charge1 > 0.F
          ? std::atan2(py2, px2)
          : std::atan2(py1, px1);

    const double signedDeltaPhi =
        wrapPhi(phiPositive - phiNegative);
    const double minimumSignedDeltaPhi =
        0.8 - 0.4 * std::min(kshortPt, 2.0);

    const bool passSelection =
        finite(momentum1) &&
        finite(momentum2) &&
        std::isfinite(opening) &&
        massKshort >= massMinimum &&
        massKshort <= massMaximum &&
        pt1 >= cut.minimumDaughterPt &&
        pt2 >= cut.minimumDaughterPt &&
        pt1 <= daughterPtMaximum &&
        pt2 <= daughterPtMaximum &&
        kshortPt <= kshortPtMaximum &&
        (optionalOpeningAngleMinimum < 0. ||
         opening >= optionalOpeningAngleMinimum) &&
        signedDeltaPhi >= minimumSignedDeltaPhi &&
        std::isfinite(dedx1) &&
        std::isfinite(dedx2) &&
        dedx1 < 400. &&
        dedx2 < 400. &&
        dira3D(totalMomentum, flight) >= cut.minimumDira3D &&
        std::abs(pairDca) <= cut.maximumPairDca &&
        std::abs(pcaZ) < cut.maximumAbsPcaZ &&
        std::abs(pca1Z - pca2Z) <
          cut.maximumAbsDeltaPcaZ &&
        decayRadius > cut.minimumDecayRadius &&
        std::abs(alpha) < cut.maximumAbsAlpha &&
        quality1 < cut.maximumQuality &&
        quality2 < cut.maximumQuality &&
        npoints1 > cut.minimumNpoints &&
        npoints2 > cut.minimumNpoints;

    if (!passSelection) continue;

    ++selectedCandidates;
    hFailureCode.Fill(0.);

    const ScaleSolution solution = solveScales(
      momentum1,
      momentum2,
      flight,
      minimumMomentumScale,
      maximumMomentumScale,
      massTolerance,
      diraXYTolerance);

    if (!solution.valid)
    {
      hFailureCode.Fill(1.);
      continue;
    }

    const int side1 = inferSide(daughter1Side, momentum1);
    const int side2 = inferSide(daughter2Side, momentum2);

    if ((side1 != 0 && side1 != 1) ||
        (side2 != 0 && side2 != 1))
    {
      hFailureCode.Fill(3.);
      continue;
    }

    ++solvedCandidates;

    const double massBefore =
        invariantMass(momentum1, momentum2);
    const double currentDiraXY =
        diraXY(totalMomentum, flight);
    const double currentDira3D =
        dira3D(totalMomentum, flight);

    hMassBefore.Fill(massBefore);
    hMassAfter.Fill(solution.massAfter);
    hDiraXYBefore.Fill(currentDiraXY);
    hDiraXYAfter.Fill(solution.diraXYAfter);
    hDira3DBefore.Fill(currentDira3D);
    hScale1.Fill(solution.scale1);
    hScale2.Fill(solution.scale2);
    hScale2VsScale1.Fill(
      solution.scale1,
      solution.scale2);

    outKshortPt = kshortPt;
    outMassBefore = massBefore;
    outMassAfter = solution.massAfter;
    outDiraXYBefore = currentDiraXY;
    outDiraXYAfter = solution.diraXYAfter;
    outDira3DBefore = currentDira3D;
    outPairDca = pairDca;
    outDecayRadius = decayRadius;
    outAlpha = alpha;

    const std::array<Vec3, 2> momenta = {
      momentum1,
      momentum2
    };

    const std::array<double, 2> scales = {
      solution.scale1,
      solution.scale2
    };

    const std::array<int, 2> sides = {
      side1,
      side2
    };

    const std::array<int, 2> charges = {
      charge1 > 0.F ? 1 : -1,
      charge2 > 0.F ? 1 : -1
    };

    const std::array<double, 2> otherPts = {
      pt2,
      pt1
    };

    const std::array<double, 2> otherScales = {
      solution.scale2,
      solution.scale1
    };

    for (int daughterIndex = 0;
         daughterIndex < 2;
         ++daughterIndex)
    {
      const Vec3& momentum = momenta[daughterIndex];

      outSide = sides[daughterIndex];
      outCharge = charges[daughterIndex];
      outDaughterIndex = daughterIndex + 1;

      outPt = transverseMomentum(momentum);
      outPhi = wrapPhi(
        std::atan2(momentum.y, momentum.x));
      outEta = pseudorapidity(momentum);
      outP = norm(momentum);
      outMomentumScale = scales[daughterIndex];

      outOtherPt = otherPts[daughterIndex];
      outOtherMomentumScale =
          otherScales[daughterIndex];

      if (!std::isfinite(outPt) ||
          !std::isfinite(outPhi) ||
          !std::isfinite(outEta) ||
          !std::isfinite(outP) ||
          !std::isfinite(outMomentumScale))
      {
        hFailureCode.Fill(2.);
        continue;
      }

      scaleTree.Fill();
      hFailureCode.Fill(4.);
      ++writtenDaughters;
    }
  }

  output->cd();

  TNamed selectionMetadata(
    "selection",
    TString::Format(
      "%s; %.3f<m(pi pi)<%.3f; pion dE/dx<400; "
      "signed DeltaPhi >= 0.8-0.4*min(KShort pT,2); "
      "optional opening angle minimum=%.6g",
      cut.name.c_str(),
      massMinimum,
      massMaximum,
      optionalOpeningAngleMinimum).Data());
  selectionMetadata.Write();

  TNamed methodMetadata(
    "method",
    "With daughter phi directions fixed, solve scale2/scale1 "
    "from exact transverse pointing and solve the remaining "
    "common scale from the PDG KShort mass.");
  methodMetadata.Write();

  output->Write();
  output->Close();

  std::cout
      << "Input entries: " << totalEntries << '\n'
      << "Selected KShort candidates: "
      << selectedCandidates << '\n'
      << "Solved candidates: "
      << solvedCandidates << '\n'
      << "Written daughter samples: "
      << writtenDaughters << '\n'
      << "Wrote: " << outputName << std::endl;
}
