// ExportK0sMomentumScaleTrainingCSV.C
//
// Export one CSV row per K0S -> pi+ pi- candidate for weakly supervised
// momentum-scale training.
//
// The CSV keeps both daughters in the same row because the training target is
// the corrected pair mass.  It also stores enough variables to reproduce the
// established K0S selections in the notebook.
//
// No momentum-scale target is assigned candidate by candidate.
//
// The recommended DNN model is trained independently for side 0 and side 1.
// Within one side, positive and negative tracks enter the same K0S mass loss.
//
// Example:
//
// root -l -b -q 'ExportK0sMomentumScaleTrainingCSV.C(
//   "/path/to/files",
//   "p_v0_*.root",
//   "output/k0s_momentum_scale_training.csv",
//   "pairTree")'

#include <TChain.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
  constexpr double kPionMass = 0.13957039;
  constexpr double kTwoPi = 2.0 * TMath::Pi();

  struct Vec3
  {
    double x = 0.;
    double y = 0.;
    double z = 0.;
  };

  Vec3 add(const Vec3& first, const Vec3& second)
  {
    return {
      first.x + second.x,
      first.y + second.y,
      first.z + second.z
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
    while (phi >= TMath::Pi())
    {
      phi -= kTwoPi;
    }

    while (phi < -TMath::Pi())
    {
      phi += kTwoPi;
    }

    return phi;
  }

  double pseudorapidity(const Vec3& momentum)
  {
    const double magnitude = norm(momentum);
    const double numerator = magnitude + momentum.z;
    const double denominator = magnitude - momentum.z;

    if (!(numerator > 0.) ||
        !(denominator > 0.))
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
      std::sqrt(
        norm2(momentum1) +
        kPionMass * kPionMass);

    const double energy2 =
      std::sqrt(
        norm2(momentum2) +
        kPionMass * kPionMass);

    const Vec3 total =
      add(momentum1, momentum2);

    const double massSquared =
      (energy1 + energy2) *
      (energy1 + energy2) -
      norm2(total);

    return massSquared > 0.
      ? std::sqrt(massSquared)
      : 0.;
  }

  double dira3D(
      const Vec3& momentum,
      const Vec3& flight)
  {
    const double denominator =
      norm(momentum) * norm(flight);

    if (!(denominator > 0.))
    {
      return -2.;
    }

    return std::clamp(
      (
        momentum.x * flight.x +
        momentum.y * flight.y +
        momentum.z * flight.z
      ) / denominator,
      -1.,
      1.);
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

    return std::clamp(
      (
        momentum.x * flight.x +
        momentum.y * flight.y
      ) /
      (
        momentumMagnitude *
        flightMagnitude
      ),
      -1.,
      1.);
  }

  int inferSide(
      const int detailedSide,
      const Vec3& momentum)
  {
    if (detailedSide == 0 ||
        detailedSide == 1)
    {
      return detailedSide;
    }

    return momentum.z >= 0.
      ? 1
      : 0;
  }

  void writeHeader(std::ofstream& output)
  {
    output
      << "entry,entry_parity,"
      << "candidate_mask,has_kshort_details,"
      << "same_side,"
      << "side1,charge1,px1,py1,pz1,pt1,phi1,eta1,p1,"
      << "side2,charge2,px2,py2,pz2,pt2,phi2,eta2,p2,"
      << "mass_kshort_branch,mass_before,kshort_pt,"
      << "pair_dca,alpha,"
      << "pca_x,pca_y,pca_z,pca1_z,pca2_z,"
      << "flight_x,flight_y,flight_z,decay_radius,"
      << "dira_xy,dira_3d,opening_angle,"
      << "dedx1,dedx2,quality1,quality2,npoints1,npoints2,"
      << "signed_delta_phi,minimum_signed_delta_phi,"
      << "pass_pion_pid,pass_signed_delta_phi,"
      << "pass_cut03,pass_cut07"
      << '\n';
  }
}

void ExportK0sMomentumScaleTrainingCSV(
    const char* inputDirectory = ".",
    const char* filePattern = "*.root",
    const char* outputCsv =
      "k0s_momentum_scale_training.csv",
    const char* treeName = "pairTree",
    const double beamX = 0.158,
    const double beamY = 0.285,
    const double beamZ = 0.0,
    const double broadMassMinimum = 0.40,
    const double broadMassMaximum = 0.60,
    const double daughterPtMaximum = 5.0,
    const double kshortPtMaximum = 3.0,
    const Long64_t maximumEntries = -1)
{
  const TString chainPattern =
    TString::Format(
      "%s/%s",
      inputDirectory,
      filePattern);

  TChain chain(treeName);

  const int numberOfFiles =
    chain.Add(chainPattern);

  if (numberOfFiles <= 0)
  {
    std::cerr
      << "ERROR: no files matched "
      << chainPattern
      << std::endl;

    return;
  }

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

  const std::vector<const char*>
    requiredBranches = {
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

  for (const char* branchName :
       requiredBranches)
  {
    if (!chain.GetBranch(branchName))
    {
      std::cerr
        << "ERROR: missing branch "
        << branchName
        << std::endl;

      return;
    }
  }

  chain.SetBranchAddress(
    "candidate_mask",
    &candidateMask);

  chain.SetBranchAddress(
    "has_kshort_daughter_details",
    &hasKshortDetails);

  chain.SetBranchAddress(
    "mass_Kshort",
    &massKshort);

  chain.SetBranchAddress(
    "alpha",
    &alpha);

  chain.SetBranchAddress(
    "pairDCA",
    &pairDca);

  chain.SetBranchAddress(
    "dedx_1",
    &dedx1);

  chain.SetBranchAddress(
    "dedx_2",
    &dedx2);

  chain.SetBranchAddress(
    "quality1",
    &quality1);

  chain.SetBranchAddress(
    "quality2",
    &quality2);

  chain.SetBranchAddress(
    "charge1",
    &charge1);

  chain.SetBranchAddress(
    "charge2",
    &charge2);

  chain.SetBranchAddress(
    "pca_x",
    &pcaX);

  chain.SetBranchAddress(
    "pca_y",
    &pcaY);

  chain.SetBranchAddress(
    "pca_z",
    &pcaZ);

  chain.SetBranchAddress(
    "pca1_z",
    &pca1Z);

  chain.SetBranchAddress(
    "pca2_z",
    &pca2Z);

  chain.SetBranchAddress(
    "px1",
    &px1);

  chain.SetBranchAddress(
    "py1",
    &py1);

  chain.SetBranchAddress(
    "pz1",
    &pz1);

  chain.SetBranchAddress(
    "px2",
    &px2);

  chain.SetBranchAddress(
    "py2",
    &py2);

  chain.SetBranchAddress(
    "pz2",
    &pz2);

  chain.SetBranchAddress(
    "npoints1",
    &npoints1);

  chain.SetBranchAddress(
    "npoints2",
    &npoints2);

  chain.SetBranchAddress(
    "daughter1_side",
    &daughter1Side);

  chain.SetBranchAddress(
    "daughter2_side",
    &daughter2Side);

  std::ofstream output(outputCsv);

  if (!output)
  {
    std::cerr
      << "ERROR: could not create "
      << outputCsv
      << std::endl;

    return;
  }

  output
    << std::setprecision(10)
    << std::scientific;

  writeHeader(output);

  const Long64_t totalEntries =
    chain.GetEntries();

  const Long64_t entriesToRun =
    maximumEntries > 0
      ? std::min(
          maximumEntries,
          totalEntries)
      : totalEntries;

  Long64_t writtenRows = 0;
  Long64_t sameSideRows = 0;

  for (Long64_t entry = 0;
       entry < entriesToRun;
       ++entry)
  {
    chain.GetEntry(entry);

    if (entry % 100000 == 0)
    {
      std::cout
        << "Processing "
        << entry
        << " / "
        << entriesToRun
        << std::endl;
    }

    if ((candidateMask & 1U) == 0U ||
        hasKshortDetails == 0 ||
        charge1 * charge2 >= 0.F)
    {
      continue;
    }

    const Vec3 momentum1{
      px1,
      py1,
      pz1
    };

    const Vec3 momentum2{
      px2,
      py2,
      pz2
    };

    if (!finite(momentum1) ||
        !finite(momentum2))
    {
      continue;
    }

    const double pt1 =
      transverseMomentum(momentum1);

    const double pt2 =
      transverseMomentum(momentum2);

    const double eta1 =
      pseudorapidity(momentum1);

    const double eta2 =
      pseudorapidity(momentum2);

    if (!std::isfinite(eta1) ||
        !std::isfinite(eta2))
    {
      continue;
    }

    const double phi1 =
      wrapPhi(
        std::atan2(
          momentum1.y,
          momentum1.x));

    const double phi2 =
      wrapPhi(
        std::atan2(
          momentum2.y,
          momentum2.x));

    const double p1 =
      norm(momentum1);

    const double p2 =
      norm(momentum2);

    const Vec3 totalMomentum =
      add(momentum1, momentum2);

    const double kshortPt =
      transverseMomentum(totalMomentum);

    const double massBefore =
      invariantMass(
        momentum1,
        momentum2);

    if (massBefore < broadMassMinimum ||
        massBefore > broadMassMaximum ||
        pt1 > daughterPtMaximum ||
        pt2 > daughterPtMaximum ||
        kshortPt > kshortPtMaximum)
    {
      continue;
    }

    const Vec3 flight{
      pcaX - beamX,
      pcaY - beamY,
      pcaZ - beamZ
    };

    const double decayRadius =
      std::hypot(
        flight.x,
        flight.y);

    const double momentumProduct =
      p1 * p2;

    const double openingAngle =
      momentumProduct > 0.
        ? std::acos(
            std::clamp(
              (
                momentum1.x *
                momentum2.x +
                momentum1.y *
                momentum2.y +
                momentum1.z *
                momentum2.z
              ) / momentumProduct,
              -1.,
              1.))
        : std::numeric_limits<double>::quiet_NaN();

    const int side1 =
      inferSide(
        daughter1Side,
        momentum1);

    const int side2 =
      inferSide(
        daughter2Side,
        momentum2);

    const int chargeSign1 =
      charge1 > 0.F
        ? 1
        : -1;

    const int chargeSign2 =
      charge2 > 0.F
        ? 1
        : -1;

    const bool sameSide =
      side1 == side2;

    const double phiPositive =
      chargeSign1 > 0
        ? phi1
        : phi2;

    const double phiNegative =
      chargeSign1 > 0
        ? phi2
        : phi1;

    const double signedDeltaPhi =
      wrapPhi(
        phiPositive -
        phiNegative);

    const double minimumSignedDeltaPhi =
      0.8 -
      0.4 *
      std::min(kshortPt, 2.0);

    const bool passPionPid =
      std::isfinite(dedx1) &&
      std::isfinite(dedx2) &&
      dedx1 < 400. &&
      dedx2 < 400.;

    const bool passSignedDeltaPhi =
      signedDeltaPhi >=
      minimumSignedDeltaPhi;

    const bool passCut03 =
      pt1 >= 0.20 &&
      pt2 >= 0.20 &&
      std::abs(pcaZ) < 15.0 &&
      std::abs(pca1Z - pca2Z) < 0.50 &&
      decayRadius > 2.0 &&
      std::abs(alpha) < 0.99 &&
      std::abs(pairDca) <= 2.00 &&
      dira3D(totalMomentum, flight) >= 0.85 &&
      quality1 < 15.0 &&
      quality2 < 15.0 &&
      npoints1 > 30 &&
      npoints2 > 30;

    const bool passCut07 =
      pt1 >= 0.20 &&
      pt2 >= 0.20 &&
      std::abs(pcaZ) < 10.0 &&
      std::abs(pca1Z - pca2Z) < 0.20 &&
      decayRadius > 2.0 &&
      std::abs(alpha) < 0.99 &&
      std::abs(pairDca) <= 0.50 &&
      dira3D(totalMomentum, flight) >= 0.95 &&
      quality1 < 10.0 &&
      quality2 < 10.0 &&
      npoints1 > 32 &&
      npoints2 > 32;

    output
      << entry << ','
      << (entry & 1LL) << ','
      << candidateMask << ','
      << static_cast<int>(hasKshortDetails) << ','
      << static_cast<int>(sameSide) << ','
      << side1 << ','
      << chargeSign1 << ','
      << momentum1.x << ','
      << momentum1.y << ','
      << momentum1.z << ','
      << pt1 << ','
      << phi1 << ','
      << eta1 << ','
      << p1 << ','
      << side2 << ','
      << chargeSign2 << ','
      << momentum2.x << ','
      << momentum2.y << ','
      << momentum2.z << ','
      << pt2 << ','
      << phi2 << ','
      << eta2 << ','
      << p2 << ','
      << massKshort << ','
      << massBefore << ','
      << kshortPt << ','
      << pairDca << ','
      << alpha << ','
      << pcaX << ','
      << pcaY << ','
      << pcaZ << ','
      << pca1Z << ','
      << pca2Z << ','
      << flight.x << ','
      << flight.y << ','
      << flight.z << ','
      << decayRadius << ','
      << diraXY(totalMomentum, flight) << ','
      << dira3D(totalMomentum, flight) << ','
      << openingAngle << ','
      << dedx1 << ','
      << dedx2 << ','
      << quality1 << ','
      << quality2 << ','
      << npoints1 << ','
      << npoints2 << ','
      << signedDeltaPhi << ','
      << minimumSignedDeltaPhi << ','
      << static_cast<int>(passPionPid) << ','
      << static_cast<int>(passSignedDeltaPhi) << ','
      << static_cast<int>(passCut03) << ','
      << static_cast<int>(passCut07)
      << '\n';

    ++writtenRows;

    if (sameSide)
    {
      ++sameSideRows;
    }
  }

  output.close();

  std::cout
    << "Input entries: "
    << totalEntries
    << '\n'
    << "Written CSV rows: "
    << writtenRows
    << '\n'
    << "Same-side rows: "
    << sameSideRows
    << '\n'
    << "Wrote: "
    << outputCsv
    << std::endl;
}
