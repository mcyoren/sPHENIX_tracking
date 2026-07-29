// ExportV0SpatialResidualTrainingCSV.C
//
// Export linked candidate and cluster CSV files from TpcV0CandidateTree output.
// K0S candidates are the training sample; a deterministic fraction of Lambda
// and anti-Lambda candidates is retained for validation. Prompt particles are
// not exported.
//
// Residual convention inherited from TpcV0CandidateTree:
//   residual = measured cluster position - fitted track position.
// Therefore the natural correction is measured - predicted_residual.
//
// This version reads pairTree only. It does not require trackTree.
//
// The daughter-detail vectors are exported only for entries where
// has_kshort_daughter_details != 0 and both daughter vector payloads are valid.
// K0S is the training channel. Lambda and anti-Lambda are retained only as
// preliminary validation samples from the subset that already has pairTree
// daughter details.
//
// Example:
// root -l -b -q 'ExportV0SpatialResidualTrainingCSV.C(
//   "/path/to/files","HITS*_V0.root","output/v0_spatial_training",
//   "pairTree",0.15,-1)'

#include <TBranch.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace
{
  constexpr double kTwoPi = 2.0 * TMath::Pi();

  struct Vec3
  {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };

  double pt(const Vec3& v) { return std::hypot(v.x, v.y); }
  double norm(const Vec3& v) { return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }

  double eta(const Vec3& v)
  {
    const double transverse = pt(v);
    return transverse > 0.0
      ? std::asinh(v.z / transverse)
      : std::numeric_limits<double>::quiet_NaN();
  }

  double wrapPhi(double value)
  {
    while (value >= TMath::Pi()) value -= kTwoPi;
    while (value < -TMath::Pi()) value += kTwoPi;
    return value;
  }

  double dira3D(const Vec3& momentum, const Vec3& flight)
  {
    const double a = norm(momentum);
    const double b = norm(flight);
    if (!(a > 0.0) || !(b > 0.0)) return -2.0;
    return std::clamp(
      (momentum.x*flight.x + momentum.y*flight.y + momentum.z*flight.z)/(a*b),
      -1.0, 1.0);
  }

  bool hasBranch(TTree* tree, const std::string& name)
  {
    return tree && tree->GetBranch(name.c_str());
  }

  template<class T>
  bool bindRequired(TTree* tree, const std::string& name, T* address)
  {
    if (!hasBranch(tree, name))
    {
      std::cerr << "ERROR: missing required branch " << name << std::endl;
      return false;
    }
    tree->SetBranchAddress(name.c_str(), address);
    return true;
  }

  template<class T>
  bool bindOptional(TTree* tree, const std::string& name, T* address)
  {
    if (!hasBranch(tree, name)) return false;
    tree->SetBranchAddress(name.c_str(), address);
    return true;
  }

  template<class T>
  bool bindVectorRequired(TTree* tree, const std::string& name, std::vector<T>** address)
  {
    if (!hasBranch(tree, name))
    {
      std::cerr << "ERROR: missing required vector branch " << name << std::endl;
      return false;
    }
    tree->SetBranchAddress(name.c_str(), address);
    return true;
  }

  template<class T>
  bool bindVectorOptional(TTree* tree, const std::string& name, std::vector<T>** address)
  {
    if (!hasBranch(tree, name)) return false;
    tree->SetBranchAddress(name.c_str(), address);
    return true;
  }

  struct DetailBranches
  {
    Int_t trackId = -1;
    Int_t charge = 0;
    Int_t side = -1;
    Int_t npoints = 0;
    UInt_t ntpcClusters = 0;
    Float_t dedx = std::numeric_limits<float>::quiet_NaN();
    Float_t px = std::numeric_limits<float>::quiet_NaN();
    Float_t py = std::numeric_limits<float>::quiet_NaN();
    Float_t pz = std::numeric_limits<float>::quiet_NaN();
    Float_t trackPt = std::numeric_limits<float>::quiet_NaN();
    Float_t trackEta = std::numeric_limits<float>::quiet_NaN();

    std::vector<unsigned int>* clusterIndex = nullptr;
    std::vector<int>* clusterSide = nullptr;
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
    std::vector<double>* residualR = nullptr;
    std::vector<double>* residualRPhi = nullptr;
    std::vector<double>* residualZ = nullptr;
    std::vector<double>* sigmaR = nullptr;
    std::vector<double>* sigmaRPhi = nullptr;
    std::vector<double>* sigmaZ = nullptr;
    std::vector<unsigned char>* recommended = nullptr;
    std::vector<double>* measurementChi2 = nullptr;
    std::vector<unsigned char>* measurementUsed = nullptr;

    bool bind(TTree* tree, const std::string& prefix, bool required, bool bindScalars = true, bool requireCartesian = true)
    {
      auto name = [&](const std::string& suffix)
      {
        return prefix.empty() ? suffix : prefix + "_" + suffix;
      };
      auto scalar = [&](const std::string& suffix, auto* address)
      {
        return required
          ? bindRequired(tree, name(suffix), address)
          : bindOptional(tree, name(suffix), address);
      };

      bool ok = true;
      if (bindScalars)
      {
        ok &= scalar("track_id", &trackId);
        ok &= scalar("charge", &charge);
        ok &= scalar("side", &side);
        ok &= scalar("npoints", &npoints);
        ok &= scalar("ntpc_clusters", &ntpcClusters);
        scalar("dedx", &dedx);
        scalar("px", &px);
        scalar("py", &py);
        scalar("pz", &pz);
        scalar("pt", &trackPt);
        scalar("eta", &trackEta);
      }

      ok &= bindVectorRequired(tree, name("cluster_index"), &clusterIndex);
      ok &= bindVectorRequired(tree, name("cluster_side"), &clusterSide);
      ok &= bindVectorRequired(tree, name("layer"), &layer);
      if (requireCartesian)
      {
        ok &= bindVectorRequired(tree, name("cluster_x"), &clusterX);
        ok &= bindVectorRequired(tree, name("cluster_y"), &clusterY);
      }
      else
      {
        bindVectorOptional(tree, name("cluster_x"), &clusterX);
        bindVectorOptional(tree, name("cluster_y"), &clusterY);
      }
      ok &= bindVectorRequired(tree, name("cluster_z"), &clusterZ);
      ok &= bindVectorRequired(tree, name("cluster_r"), &clusterR);
      ok &= bindVectorRequired(tree, name("cluster_phi"), &clusterPhi);
      ok &= bindVectorRequired(tree, name("residual_r"), &residualR);
      ok &= bindVectorRequired(tree, name("residual_rphi"), &residualRPhi);
      ok &= bindVectorRequired(tree, name("residual_z"), &residualZ);

      bindVectorOptional(tree, name("fit_x"), &fitX);
      bindVectorOptional(tree, name("fit_y"), &fitY);
      bindVectorOptional(tree, name("fit_z"), &fitZ);
      bindVectorOptional(tree, name("fit_px"), &fitPx);
      bindVectorOptional(tree, name("fit_py"), &fitPy);
      bindVectorOptional(tree, name("fit_pz"), &fitPz);
      bindVectorOptional(tree, name("assigned_sigma_r"), &sigmaR);
      bindVectorOptional(tree, name("assigned_sigma_rphi"), &sigmaRPhi);
      bindVectorOptional(tree, name("assigned_sigma_z"), &sigmaZ);
      bindVectorOptional(tree, name("recommended_for_reference_fit"), &recommended);
      bindVectorOptional(tree, name("kalman_measurement_chi2"), &measurementChi2);
      bindVectorOptional(tree, name("kalman_measurement_used"), &measurementUsed);
      return ok;
    }
  };

  struct TrackDetail
  {
    int trackId = -1;
    int charge = 0;
    int side = -1;
    int npoints = 0;
    unsigned int ntpcClusters = 0;
    double dedx = std::numeric_limits<double>::quiet_NaN();
    double px = std::numeric_limits<double>::quiet_NaN();
    double py = std::numeric_limits<double>::quiet_NaN();
    double pz = std::numeric_limits<double>::quiet_NaN();
    double trackPt = std::numeric_limits<double>::quiet_NaN();
    double trackEta = std::numeric_limits<double>::quiet_NaN();

    std::vector<unsigned int> clusterIndex;
    std::vector<int> clusterSide;
    std::vector<unsigned int> layer;
    std::vector<double> clusterX, clusterY, clusterZ, clusterR, clusterPhi;
    std::vector<double> fitX, fitY, fitZ, fitPx, fitPy, fitPz;
    std::vector<double> residualR, residualRPhi, residualZ;
    std::vector<double> sigmaR, sigmaRPhi, sigmaZ;
    std::vector<unsigned char> recommended;
    std::vector<double> measurementChi2;
    std::vector<unsigned char> measurementUsed;

    bool valid() const
    {
      const std::size_t n = clusterR.size();
      return n > 0 && clusterPhi.size() == n && clusterZ.size() == n &&
             residualR.size() == n && residualRPhi.size() == n && residualZ.size() == n;
    }
  };

  template<class T>
  std::vector<T> copyVector(const std::vector<T>* input)
  {
    return input ? *input : std::vector<T>{};
  }

  TrackDetail copyDetail(const DetailBranches& input)
  {
    TrackDetail out;
    out.trackId = input.trackId;
    out.charge = input.charge;
    out.side = input.side;
    out.npoints = input.npoints;
    out.ntpcClusters = input.ntpcClusters;
    out.dedx = input.dedx;
    out.px = input.px;
    out.py = input.py;
    out.pz = input.pz;
    out.trackPt = input.trackPt;
    out.trackEta = input.trackEta;
    out.clusterIndex = copyVector(input.clusterIndex);
    out.clusterSide = copyVector(input.clusterSide);
    out.layer = copyVector(input.layer);
    out.clusterX = copyVector(input.clusterX);
    out.clusterY = copyVector(input.clusterY);
    out.clusterZ = copyVector(input.clusterZ);
    out.clusterR = copyVector(input.clusterR);
    out.clusterPhi = copyVector(input.clusterPhi);
    out.fitX = copyVector(input.fitX);
    out.fitY = copyVector(input.fitY);
    out.fitZ = copyVector(input.fitZ);
    out.fitPx = copyVector(input.fitPx);
    out.fitPy = copyVector(input.fitPy);
    out.fitPz = copyVector(input.fitPz);
    out.residualR = copyVector(input.residualR);
    out.residualRPhi = copyVector(input.residualRPhi);
    out.residualZ = copyVector(input.residualZ);
    out.sigmaR = copyVector(input.sigmaR);
    out.sigmaRPhi = copyVector(input.sigmaRPhi);
    out.sigmaZ = copyVector(input.sigmaZ);
    out.recommended = copyVector(input.recommended);
    out.measurementChi2 = copyVector(input.measurementChi2);
    out.measurementUsed = copyVector(input.measurementUsed);
    return out;
  }

  template<class T>
  T atOr(const std::vector<T>& values, std::size_t index, T fallback)
  {
    return index < values.size() ? values[index] : fallback;
  }

  std::uint64_t hash64(int run, int event, Long64_t entry, int channel)
  {
    std::uint64_t x = static_cast<std::uint32_t>(run);
    x ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(event)) << 32;
    x ^= static_cast<std::uint64_t>(entry) * 0x9E3779B97F4A7C15ULL;
    x ^= static_cast<std::uint64_t>(channel + 1) * 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27; x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;
    return x;
  }

  double uniformHash(int run, int event, Long64_t entry, int channel)
  {
    return static_cast<double>(hash64(run, event, entry, channel) >> 11) /
           9007199254740992.0;
  }

  int particleHypothesis(int channel, int charge)
  {
    if (channel == 0) return charge > 0 ? 211 : -211;
    if (channel == 1) return charge > 0 ? 2212 : -211;
    if (channel == 2) return charge < 0 ? -2212 : 211;
    return 0;
  }

  void writeCandidateHeader(std::ofstream& out)
  {
    out << "candidate_id,file_index,pair_entry,run,evt,event_parity,channel,candidate_mask,"
        << "track_id1,track_id2,charge1,charge2,side1,side2,same_side,"
        << "px1,py1,pz1,pt1,phi1,eta1,px2,py2,pz2,pt2,phi2,eta2,"
        << "pca_x,pca_y,pca_z,pca1_z,pca2_z,pair_dca,alpha,dira,decay_radius,"
        << "dedx1,dedx2,quality1,quality2,npoints1,npoints2,ntpc_clusters1,ntpc_clusters2,"
        << "mass_kshort,mass_lambda,mass_antilambda,selected_mass,pair_pt,"
        << "pass_cut03,pass_pion_pid,pass_signed_delta_phi,"
        << "number_of_clusters1,number_of_clusters2\n";
  }

  void writeClusterHeader(std::ofstream& out)
  {
    out << "candidate_id,channel,daughter,particle_hypothesis,run,evt,event_parity,"
        << "track_id,charge,side,track_px,track_py,track_pz,track_pt,track_phi,track_eta,"
        << "track_npoints,track_ntpc_clusters,cluster_ordinal,cluster_index,layer,cluster_side,"
        << "cluster_x,cluster_y,cluster_z,cluster_r,cluster_phi,"
        << "fit_x,fit_y,fit_z,fit_px,fit_py,fit_pz,"
        << "residual_r,residual_rphi,residual_z,"
        << "assigned_sigma_r,assigned_sigma_rphi,assigned_sigma_z,"
        << "recommended_for_reference_fit,kalman_measurement_chi2,kalman_measurement_used\n";
  }

  void writeClusters(
      std::ofstream& out,
      Long64_t candidateId,
      int channel,
      int daughter,
      int run,
      int event,
      int parity,
      const TrackDetail& track)
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double trackPt = std::isfinite(track.trackPt)
      ? track.trackPt : std::hypot(track.px, track.py);
    const double trackPhi = std::atan2(track.py, track.px);
    const double trackEta = std::isfinite(track.trackEta)
      ? track.trackEta : eta({track.px, track.py, track.pz});

    for (std::size_t i = 0; i < track.clusterR.size(); ++i)
    {
      const double r = atOr(track.clusterR, i, nan);
      const double phi = atOr(track.clusterPhi, i, nan);
      const double x = atOr(track.clusterX, i,
        std::isfinite(r) && std::isfinite(phi) ? r*std::cos(phi) : nan);
      const double y = atOr(track.clusterY, i,
        std::isfinite(r) && std::isfinite(phi) ? r*std::sin(phi) : nan);

      out << candidateId << ',' << channel << ',' << daughter << ','
          << particleHypothesis(channel, track.charge) << ','
          << run << ',' << event << ',' << parity << ','
          << track.trackId << ',' << track.charge << ',' << track.side << ','
          << track.px << ',' << track.py << ',' << track.pz << ','
          << trackPt << ',' << trackPhi << ',' << trackEta << ','
          << track.npoints << ',' << track.ntpcClusters << ','
          << i << ',' << atOr(track.clusterIndex, i, 0U) << ','
          << atOr(track.layer, i, 0U) << ',' << atOr(track.clusterSide, i, track.side) << ','
          << x << ',' << y << ',' << atOr(track.clusterZ, i, nan) << ','
          << r << ',' << phi << ','
          << atOr(track.fitX, i, nan) << ',' << atOr(track.fitY, i, nan) << ','
          << atOr(track.fitZ, i, nan) << ',' << atOr(track.fitPx, i, nan) << ','
          << atOr(track.fitPy, i, nan) << ',' << atOr(track.fitPz, i, nan) << ','
          << atOr(track.residualR, i, nan) << ','
          << atOr(track.residualRPhi, i, nan) << ','
          << atOr(track.residualZ, i, nan) << ','
          << atOr(track.sigmaR, i, nan) << ','
          << atOr(track.sigmaRPhi, i, nan) << ','
          << atOr(track.sigmaZ, i, nan) << ','
          << static_cast<int>(atOr(track.recommended, i, static_cast<unsigned char>(1))) << ','
          << atOr(track.measurementChi2, i, nan) << ','
          << static_cast<int>(atOr(track.measurementUsed, i, static_cast<unsigned char>(1)))
          << '\n';
    }
  }
}

// NOTE:
// Lambda and anti-Lambda output from an existing file is restricted to the
// subset for which pairTree already contains daughter-detail vectors. This
// sample is suitable for a preliminary cross-check, but it is not guaranteed
// to be an unbiased or complete Lambda validation sample until the producer
// fills daughter details for all selected V0 channels.
void ExportV0SpatialResidualTrainingCSV(
    const char* inputDirectory = ".",
    const char* filePattern = "*.root",
    const char* outputPrefix = "v0_spatial_training",
    const char* pairTreeName = "pairTree",
    const double lambdaKeepFraction = 1.0,
    const Long64_t maximumPairEntries = -1,
    const double beamX = 0.158,
    const double beamY = 0.285,
    const double beamZ = 0.0,
    const double kshortMassMinimum = 0.40,
    const double kshortMassMaximum = 0.60,
    const double lambdaMassMinimum = 1.05,
    const double lambdaMassMaximum = 1.25)
{
  TChain discovery(pairTreeName);
  const TString pattern = TString::Format("%s/%s", inputDirectory, filePattern);
  const int nfiles = discovery.Add(pattern);
  if (nfiles <= 0)
  {
    std::cerr << "ERROR: no files matched " << pattern << std::endl;
    return;
  }

  std::ofstream candidateOut(std::string(outputPrefix) + "_candidates.csv");
  std::ofstream clusterOut(std::string(outputPrefix) + "_clusters.csv");
  if (!candidateOut || !clusterOut)
  {
    std::cerr << "ERROR: cannot create CSV outputs for prefix " << outputPrefix << std::endl;
    return;
  }
  candidateOut << std::setprecision(10);
  clusterOut << std::setprecision(10);
  writeCandidateHeader(candidateOut);
  writeClusterHeader(clusterOut);

  Long64_t candidateId = 0;
  Long64_t processed = 0;
  Long64_t nKshort = 0;
  Long64_t nLambda = 0;
  Long64_t nAntiLambda = 0;
  Long64_t nSelectedKshortWithoutDetails = 0;
  Long64_t nSelectedLambdaWithoutDetails = 0;
  Long64_t nSelectedAntiLambdaWithoutDetails = 0;
  Long64_t nInvalidDaughterVectors = 0;

  TObjArray* files = discovery.GetListOfFiles();
  for (int fileIndex = 0; fileIndex < files->GetEntries(); ++fileIndex)
  {
    if (maximumPairEntries >= 0 && processed >= maximumPairEntries) break;

    auto* element = dynamic_cast<TChainElement*>(files->At(fileIndex));
    if (!element) continue;

    const std::string fileName = element->GetTitle();
    std::unique_ptr<TFile> input(TFile::Open(fileName.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      std::cerr << "WARNING: cannot open " << fileName << std::endl;
      continue;
    }

    TTree* pairTree = dynamic_cast<TTree*>(input->Get(pairTreeName));
    if (!pairTree)
    {
      std::cerr << "WARNING: missing " << pairTreeName << " in " << fileName << std::endl;
      continue;
    }

    Int_t run = 0, evt = 0;
    UInt_t candidateMask = 0;
    Int_t hasDetails = 0;
    Float_t px1 = 0, py1 = 0, pz1 = 0, px2 = 0, py2 = 0, pz2 = 0;
    Float_t pcaX = 0, pcaY = 0, pcaZ = 0, pca1Z = 0, pca2Z = 0;
    Float_t pairDca = 0, alpha = 0, charge1 = 0, charge2 = 0;
    Float_t dedx1 = 0, dedx2 = 0, quality1 = 0, quality2 = 0;
    Int_t trackId1 = -1, trackId2 = -1;
    Short_t npoints1 = 0, npoints2 = 0;
    Float_t massKshort = 0, massLambda = 0, massAntiLambda = 0;

    bool ok = true;
    ok &= bindRequired(pairTree, "run", &run);
    ok &= bindRequired(pairTree, "evt", &evt);
    ok &= bindRequired(pairTree, "candidate_mask", &candidateMask);
    ok &= bindRequired(pairTree, "has_kshort_daughter_details", &hasDetails);
    ok &= bindRequired(pairTree, "px1", &px1);
    ok &= bindRequired(pairTree, "py1", &py1);
    ok &= bindRequired(pairTree, "pz1", &pz1);
    ok &= bindRequired(pairTree, "px2", &px2);
    ok &= bindRequired(pairTree, "py2", &py2);
    ok &= bindRequired(pairTree, "pz2", &pz2);
    ok &= bindRequired(pairTree, "pca_x", &pcaX);
    ok &= bindRequired(pairTree, "pca_y", &pcaY);
    ok &= bindRequired(pairTree, "pca_z", &pcaZ);
    ok &= bindRequired(pairTree, "pca1_z", &pca1Z);
    ok &= bindRequired(pairTree, "pca2_z", &pca2Z);
    ok &= bindRequired(pairTree, "pairDCA", &pairDca);
    ok &= bindRequired(pairTree, "alpha", &alpha);
    ok &= bindRequired(pairTree, "charge1", &charge1);
    ok &= bindRequired(pairTree, "charge2", &charge2);
    ok &= bindRequired(pairTree, "dedx_1", &dedx1);
    ok &= bindRequired(pairTree, "dedx_2", &dedx2);
    ok &= bindRequired(pairTree, "quality1", &quality1);
    ok &= bindRequired(pairTree, "quality2", &quality2);
    ok &= bindRequired(pairTree, "track_id1", &trackId1);
    ok &= bindRequired(pairTree, "track_id2", &trackId2);
    ok &= bindRequired(pairTree, "npoints1", &npoints1);
    ok &= bindRequired(pairTree, "npoints2", &npoints2);
    ok &= bindRequired(pairTree, "mass_Kshort", &massKshort);
    ok &= bindRequired(pairTree, "mass_Lambda", &massLambda);
    ok &= bindRequired(pairTree, "mass_AntiLambda", &massAntiLambda);

    DetailBranches pairDetail1;
    DetailBranches pairDetail2;

    const bool pairDetailsAvailable =
      hasBranch(pairTree, "daughter1_cluster_r") &&
      hasBranch(pairTree, "daughter2_cluster_r") &&
      hasBranch(pairTree, "daughter1_residual_r") &&
      hasBranch(pairTree, "daughter2_residual_r") &&
      hasBranch(pairTree, "daughter1_residual_rphi") &&
      hasBranch(pairTree, "daughter2_residual_rphi") &&
      hasBranch(pairTree, "daughter1_residual_z") &&
      hasBranch(pairTree, "daughter2_residual_z");

    if (!pairDetailsAvailable)
    {
      std::cerr
        << "ERROR: pairTree in "
        << fileName
        << " does not contain the required daughter-detail vectors"
        << std::endl;

      continue;
    }

    ok &= pairDetail1.bind(
      pairTree,
      "daughter1",
      true);

    ok &= pairDetail2.bind(
      pairTree,
      "daughter2",
      true);

    if (!ok)
    {
      std::cerr
        << "ERROR: branch binding failed in "
        << fileName
        << std::endl;

      continue;
    }

    const Long64_t entries = pairTree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry)
    {
      if (maximumPairEntries >= 0 && processed >= maximumPairEntries) break;
      pairTree->GetEntry(entry);
      ++processed;
      if (processed % 100000 == 0)
        std::cout << "Processed pair entries: " << processed << std::endl;

      const Vec3 p1{px1, py1, pz1};
      const Vec3 p2{px2, py2, pz2};
      const Vec3 pair{px1 + px2, py1 + py2, pz1 + pz2};
      const Vec3 flight{pcaX - beamX, pcaY - beamY, pcaZ - beamZ};
      const double pt1 = pt(p1), pt2 = pt(p2);
      const double eta1 = eta(p1), eta2 = eta(p2);
      const double phi1 = std::atan2(py1, px1), phi2 = std::atan2(py2, px2);
      const double pairPt = pt(pair);
      const double decayRadius = std::hypot(flight.x, flight.y);
      const double dira = dira3D(pair, flight);
      const bool unlike = charge1*charge2 < 0.0;
      const double deltaPcaZ = std::abs(pca1Z - pca2Z);
      const bool passCut03 = unlike && std::abs(pcaZ) < 15.0 && deltaPcaZ < 0.50 &&
        std::min(pt1, pt2) > 0.20 && decayRadius > 2.0 && std::abs(alpha) < 0.99 &&
        std::abs(pairDca) < 2.0 && dira > 0.85 &&
        std::max<double>(quality1, quality2) < 15.0 &&
        std::min<int>(npoints1, npoints2) > 30;
      const bool passPionPid = dedx1 < 400.0 && dedx2 < 400.0;

      const double phiPositive = charge1 > 0.0 ? phi1 : phi2;
      const double phiNegative = charge1 > 0.0 ? phi2 : phi1;
      const double signedDeltaPhi = wrapPhi(phiPositive - phiNegative);
      const bool passSignedDeltaPhi = signedDeltaPhi >= 0.8 - 0.4*std::min(pairPt, 2.0);

      const bool positiveHigher = charge1 > 0.0 ? pt1 > pt2 : pt2 > pt1;
      const bool negativeHigher = charge1 > 0.0 ? pt2 > pt1 : pt1 > pt2;

      const bool selectKshort = (candidateMask & 1U) && passCut03 && passPionPid &&
        passSignedDeltaPhi && std::isfinite(massKshort) &&
        massKshort >= kshortMassMinimum && massKshort <= kshortMassMaximum;
      const bool selectLambda = (candidateMask & 2U) && passCut03 && positiveHigher &&
        std::isfinite(massLambda) && massLambda >= lambdaMassMinimum &&
        massLambda <= lambdaMassMaximum && uniformHash(run, evt, entry, 1) < lambdaKeepFraction;
      const bool selectAntiLambda = (candidateMask & 4U) && passCut03 && negativeHigher &&
        std::isfinite(massAntiLambda) && massAntiLambda >= lambdaMassMinimum &&
        massAntiLambda <= lambdaMassMaximum && uniformHash(run, evt, entry, 2) < lambdaKeepFraction;

      if (!selectKshort && !selectLambda && !selectAntiLambda) continue;

      // This exporter deliberately has no trackTree fallback.
      if (hasDetails == 0)
      {
        if (selectKshort)
        {
          ++nSelectedKshortWithoutDetails;
        }

        if (selectLambda)
        {
          ++nSelectedLambdaWithoutDetails;
        }

        if (selectAntiLambda)
        {
          ++nSelectedAntiLambdaWithoutDetails;
        }

        continue;
      }

      TrackDetail first =
        copyDetail(pairDetail1);

      TrackDetail second =
        copyDetail(pairDetail2);

      if (!first.valid() ||
          !second.valid())
      {
        ++nInvalidDaughterVectors;
        continue;
      }

      // Pair-level values are the daughter momenta at the secondary PCA.
      // Daughter vectors provide clusters, fitted states, and residuals.
      first.trackId = trackId1;
      first.charge = static_cast<int>(std::lround(charge1));
      first.npoints = npoints1;
      first.ntpcClusters = static_cast<unsigned int>(first.clusterR.size());
      first.dedx = dedx1;
      first.px = px1;
      first.py = py1;
      first.pz = pz1;
      first.trackPt = pt1;
      first.trackEta = eta1;

      second.trackId = trackId2;
      second.charge = static_cast<int>(std::lround(charge2));
      second.npoints = npoints2;
      second.ntpcClusters = static_cast<unsigned int>(second.clusterR.size());
      second.dedx = dedx2;
      second.px = px2;
      second.py = py2;
      second.pz = pz2;
      second.trackPt = pt2;
      second.trackEta = eta2;

      const int side1 = (first.side == 0 || first.side == 1) ? first.side : (eta1 >= 0.0 ? 1 : 0);
      const int side2 = (second.side == 0 || second.side == 1) ? second.side : (eta2 >= 0.0 ? 1 : 0);
      first.side = side1;
      second.side = side2;
      const int parity = static_cast<int>(hash64(run, evt, 0, 0) & 1ULL);

      const int channels[3] = {selectKshort ? 0 : -1, selectLambda ? 1 : -1,
                               selectAntiLambda ? 2 : -1};
      for (int channel : channels)
      {
        if (channel < 0) continue;
        const double selectedMass = channel == 0 ? massKshort :
          (channel == 1 ? massLambda : massAntiLambda);

        candidateOut
          << candidateId << ',' << fileIndex << ',' << entry << ',' << run << ',' << evt << ','
          << parity << ',' << channel << ',' << candidateMask << ',' << trackId1 << ',' << trackId2 << ','
          << charge1 << ',' << charge2 << ',' << side1 << ',' << side2 << ',' << int(side1 == side2) << ','
          << px1 << ',' << py1 << ',' << pz1 << ',' << pt1 << ',' << phi1 << ',' << eta1 << ','
          << px2 << ',' << py2 << ',' << pz2 << ',' << pt2 << ',' << phi2 << ',' << eta2 << ','
          << pcaX << ',' << pcaY << ',' << pcaZ << ',' << pca1Z << ',' << pca2Z << ','
          << pairDca << ',' << alpha << ',' << dira << ',' << decayRadius << ','
          << dedx1 << ',' << dedx2 << ',' << quality1 << ',' << quality2 << ','
          << npoints1 << ',' << npoints2 << ',' << first.ntpcClusters << ',' << second.ntpcClusters << ','
          << massKshort << ',' << massLambda << ',' << massAntiLambda << ',' << selectedMass << ',' << pairPt << ','
          << int(passCut03) << ',' << int(passPionPid) << ',' << int(passSignedDeltaPhi) << ','
          << first.clusterR.size() << ',' << second.clusterR.size() << '\n';

        writeClusters(clusterOut, candidateId, channel, 1, run, evt, parity, first);
        writeClusters(clusterOut, candidateId, channel, 2, run, evt, parity, second);
        if (channel == 0) ++nKshort;
        else if (channel == 1) ++nLambda;
        else ++nAntiLambda;
        ++candidateId;
      }
    }
  }

  candidateOut.close();
  clusterOut.close();

  std::cout << "Wrote " << outputPrefix << "_candidates.csv and "
            << outputPrefix << "_clusters.csv" << std::endl;
  std::cout
    << "K0S: "
    << nKshort
    << ", Lambda: "
    << nLambda
    << ", anti-Lambda: "
    << nAntiLambda
    << std::endl;

  std::cout
    << "Selected entries without pairTree daughter details:"
    << " K0S="
    << nSelectedKshortWithoutDetails
    << ", Lambda="
    << nSelectedLambdaWithoutDetails
    << ", anti-Lambda="
    << nSelectedAntiLambdaWithoutDetails
    << std::endl;

  std::cout
    << "Entries with has_kshort_daughter_details=1 "
    << "but invalid or inconsistent vectors: "
    << nInvalidDaughterVectors
    << std::endl;
}
