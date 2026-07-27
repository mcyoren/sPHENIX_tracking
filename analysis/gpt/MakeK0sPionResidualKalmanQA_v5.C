// MakeK0sPionResidualKalmanQA.C
//
// One-pass QA of pion daughter tracks from K0S candidates in pairTree.
//
// Default signal region:
//   0.47 <= mass_Kshort <= 0.53 GeV/c^2
//
// Default outside-mass sideband:
//   0.40 <= mass_Kshort < 0.47  OR  0.53 < mass_Kshort <= 0.60 GeV/c^2
//
// Two V0 selections are produced:
//   cut03_baseline
//   cut07_pairDCA_5mm
//
// Important:
//   * pT > 0.20 GeV/c is used for both cuts.
//   * No primary-track DCA cut is applied by default because K0S daughters
//     are displaced tracks. DCAxy and DCAz are plotted as QA variables.
//   * Cluster residuals are rejected independently only when |residual|
//     exceeds maxAbsResidualCm (default 2 cm).
//   * daughter*_kalman_measurement_chi2 and daughter*_kalman_measurement_used
//     are filled versus Kalman measurement index. The current tree does not
//     store the original-index mapping needed to associate them unambiguously
//     with detector layer.
//
// Example:
// root -l -b -q \
// 'MakeK0sPionResidualKalmanQA_v2.C("/path/to/files","*V0*.root","output","k0s_pion_qa_v2.root")'

#include <TBranch.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace
{
  constexpr unsigned int kCandidateKShort = 1U << 0;

  struct V0Cut
  {
    std::string name;
    std::string description;
    double maxAbsPcaZ;
    double maxAbsDeltaPcaZ;
    double minDaughterPt;
    double minDecayRadius;
    double maxAbsAlpha;
    double maxAbsPairDCA;
    double minDIRA;
    double maxQuality;
    int minNpoints;
  };

  struct DaughterBranches
  {
    Int_t track_id = 0;
    Int_t charge = 0;
    Int_t side = -1;
    Int_t npoints = 0;
    UInt_t ntpc_clusters = 0;

    Float_t dedx = 0.F;
    Float_t px = 0.F;
    Float_t py = 0.F;
    Float_t pz = 0.F;
    Float_t pt = 0.F;
    Float_t eta = 0.F;
    Float_t fit_chi2 = 0.F;
    Int_t fit_ndf = 0;
    Float_t fit_chi2_ndf = 0.F;

    std::vector<unsigned int>* cluster_index = nullptr;
    std::vector<int>* cluster_side = nullptr;
    std::vector<unsigned int>* layer = nullptr;

    std::vector<double>* cluster_x = nullptr;
    std::vector<double>* cluster_y = nullptr;
    std::vector<double>* cluster_z = nullptr;
    std::vector<double>* cluster_r = nullptr;
    std::vector<double>* cluster_phi = nullptr;

    std::vector<double>* fit_x = nullptr;
    std::vector<double>* fit_y = nullptr;
    std::vector<double>* fit_z = nullptr;
    std::vector<double>* fit_px = nullptr;
    std::vector<double>* fit_py = nullptr;
    std::vector<double>* fit_pz = nullptr;

    std::vector<double>* residual_r = nullptr;
    std::vector<double>* residual_rphi = nullptr;
    std::vector<double>* residual_z = nullptr;

    std::vector<double>* assigned_sigma_r = nullptr;
    std::vector<double>* assigned_sigma_rphi = nullptr;
    std::vector<double>* assigned_sigma_z = nullptr;

    std::vector<unsigned char>* recommended_for_reference_fit = nullptr;
    std::vector<double>* kalman_measurement_chi2 = nullptr;
    std::vector<unsigned char>* kalman_measurement_used = nullptr;
  };

  struct PairHistSet
  {
    TH1D* h_mass = nullptr;
    TH1D* h_v0_pt = nullptr;
    TH1D* h_pca_z = nullptr;
    TH1D* h_abs_delta_pca_z = nullptr;
    TH1D* h_pair_dca = nullptr;
    TH1D* h_decay_radius = nullptr;
    TH1D* h_dira = nullptr;
    TH1D* h_alpha = nullptr;
    TH1D* h_qt = nullptr;

    TH2D* h_mass_vs_v0_pt = nullptr;
    TH2D* h_pair_dca_vs_delta_pca_z = nullptr;
    TH2D* h_daughter_pt_correlation = nullptr;
    TH2D* h_dca_xy_correlation = nullptr;
    TH2D* h_dca_z_correlation = nullptr;
  };

  struct DaughterHistSet
  {
    // Track-level QA.
    TH1D* h_pt = nullptr;
    TH1D* h_eta = nullptr;
    TH1D* h_phi = nullptr;
    TH1D* h_p = nullptr;
    TH1D* h_q_over_pt = nullptr;
    TH1D* h_npoints = nullptr;
    TH1D* h_ntpc_clusters = nullptr;
    TH1D* h_fit_chi2_ndf = nullptr;
    TH1D* h_dca_xy_to_primary_vertex = nullptr;
    TH1D* h_dca_z_to_primary_vertex = nullptr;
    TH1D* h_number_of_detail_points = nullptr;

    TH2D* h_pt_vs_eta = nullptr;
    TH2D* h_dca_xy_vs_phi = nullptr;
    TH2D* h_dca_z_vs_phi = nullptr;
    TH2D* h_dedx_vs_signed_p = nullptr;
    TH2D* h_npoints_vs_ntpc_clusters = nullptr;
    TH2D* h_fit_chi2_ndf_vs_pt = nullptr;

    // Track-level phase-space QA with pT as the third axis.
    TH3D* h_primary_dca_xy_vs_phi_vs_pt = nullptr;
    TH3D* h_primary_dca_z_vs_eta_vs_pt = nullptr;
    TH2D* h_secondary_dca_xy_vs_q_over_pt = nullptr;
    TH2D* h_secondary_dca_z_vs_q_over_pt = nullptr;
    TH3D* h_secondary_dca_xy_vs_phi_vs_pt = nullptr;
    TH3D* h_secondary_dca_xy_vs_eta_vs_pt = nullptr;
    TH3D* h_secondary_dca_z_vs_phi_vs_pt = nullptr;
    TH3D* h_secondary_dca_z_vs_eta_vs_pt = nullptr;
    TProfile2D* p_secondary_dca_xy_phi_eta_pt0p2 = nullptr;
    TProfile2D* p_secondary_dca_xy_phi_eta_pt1 = nullptr;
    TProfile2D* p_secondary_dca_xy_phi_eta_pt2 = nullptr;
    TProfile2D* p_secondary_dca_z_phi_eta_pt0p2 = nullptr;
    TProfile2D* p_secondary_dca_z_phi_eta_pt1 = nullptr;
    TProfile2D* p_secondary_dca_z_phi_eta_pt2 = nullptr;
    TH3D* h_quality_vs_phi_vs_pt = nullptr;
    TH3D* h_quality_vs_eta_vs_pt = nullptr;
    TH3D* h_phi_vs_eta_vs_pt = nullptr;

    // Cluster residuals and assigned measurement uncertainties.
    TH1D* h_residual_r = nullptr;
    TH1D* h_residual_rphi = nullptr;
    TH1D* h_residual_z = nullptr;
    TH1D* h_pull_r = nullptr;
    TH1D* h_pull_rphi = nullptr;
    TH1D* h_pull_z = nullptr;

    TH2D* h_residual_r_vs_layer = nullptr;
    TH2D* h_residual_rphi_vs_layer = nullptr;
    TH2D* h_residual_z_vs_layer = nullptr;
    TH2D* h_pull_r_vs_layer = nullptr;
    TH2D* h_pull_rphi_vs_layer = nullptr;
    TH2D* h_pull_z_vs_layer = nullptr;

    TH2D* h_residual_r_vs_cluster_r = nullptr;
    TH2D* h_residual_rphi_vs_cluster_r = nullptr;
    TH2D* h_residual_z_vs_cluster_r = nullptr;
    TH2D* h_residual_rphi_vs_cluster_phi = nullptr;
    TH2D* h_residual_z_vs_cluster_phi = nullptr;
    TH2D* h_residual_z_vs_cluster_z = nullptr;
    TH3D* h_residual_r_vs_cluster_r_vs_pt = nullptr;
    TH3D* h_residual_rphi_vs_cluster_r_vs_pt = nullptr;
    TH3D* h_residual_z_vs_cluster_r_vs_pt = nullptr;
    TH3D* h_residual_z_vs_cluster_z_vs_pt = nullptr;

    // Residual phase-space histograms with daughter pT as the third axis.
    TH3D* h_residual_r_vs_phi_vs_pt = nullptr;
    TH3D* h_residual_rphi_vs_phi_vs_pt = nullptr;
    TH3D* h_residual_z_vs_phi_vs_pt = nullptr;
    TH3D* h_residual_r_vs_eta_vs_pt = nullptr;
    TH3D* h_residual_rphi_vs_eta_vs_pt = nullptr;
    TH3D* h_residual_z_vs_eta_vs_pt = nullptr;
    TH3D* h_pull_rphi_vs_phi_vs_pt = nullptr;
    TH3D* h_pull_z_vs_eta_vs_pt = nullptr;
    TH3D* h_residual_rphi_vs_layer_vs_pt = nullptr;
    TH3D* h_residual_z_vs_layer_vs_pt = nullptr;

    TH2D* h_assigned_sigma_r_vs_layer = nullptr;
    TH2D* h_assigned_sigma_rphi_vs_layer = nullptr;
    TH2D* h_assigned_sigma_z_vs_layer = nullptr;
    TH2D* h_cluster_occupancy_layer_phi = nullptr;
    TH2D* h_cluster_occupancy_cluster_z_layer = nullptr;

    TProfile* p_mean_residual_r_vs_layer = nullptr;
    TProfile* p_mean_residual_rphi_vs_layer = nullptr;
    TProfile* p_mean_residual_z_vs_layer = nullptr;
    TProfile* p_recommended_fraction_vs_layer = nullptr;

    TProfile2D* p_mean_residual_rphi_layer_phi = nullptr;
    TProfile2D* p_mean_residual_z_layer_phi = nullptr;
    TProfile2D* p_mean_residual_z_cluster_z_layer = nullptr;
    TProfile2D* p_residual_rphi_phi_eta_pt0p2 = nullptr;
    TProfile2D* p_residual_rphi_phi_eta_pt1 = nullptr;
    TProfile2D* p_residual_rphi_phi_eta_pt2 = nullptr;
    TProfile2D* p_residual_z_phi_eta_pt0p2 = nullptr;
    TProfile2D* p_residual_z_phi_eta_pt1 = nullptr;
    TProfile2D* p_residual_z_phi_eta_pt2 = nullptr;
    TProfile2D* p_mean_residual_rphi_layer_phi_reference = nullptr;
    TProfile2D* p_mean_residual_z_layer_phi_reference = nullptr;

    // Fitted momentum along the daughter trajectory.
    TH2D* h_fit_pt_vs_layer = nullptr;
    TH2D* h_fit_q_over_pt_vs_layer = nullptr;
    TH2D* h_fit_pt_over_track_pt_vs_layer = nullptr;
    TH2D* h_fit_pz_vs_layer = nullptr;

    // Kalman measurement QA. These use measurement index, not detector layer.
    TH1D* h_kalman_measurement_chi2 = nullptr;
    TH1D* h_kalman_measurement_chi2_used = nullptr;
    TH1D* h_kalman_measurement_chi2_rejected = nullptr;
    TH1D* h_kalman_used_fraction = nullptr;
    TH1D* h_kalman_number_of_measurements = nullptr;
    TH1D* h_kalman_size_minus_detail_points = nullptr;
    TH2D* h_kalman_measurement_chi2_vs_index = nullptr;
    TH2D* h_kalman_measurement_used_vs_index = nullptr;
    TH2D* h_kalman_measurement_chi2_vs_pt = nullptr;
    TH3D* h_kalman_measurement_chi2_vs_index_vs_pt = nullptr;
    TH3D* h_kalman_measurement_used_vs_index_vs_pt = nullptr;
  };

  struct SampleHistSet
  {
    PairHistSet pair;
    std::map<std::string, DaughterHistSet> daughters;
  };

  bool branchExists(TChain& chain, const std::string& name)
  {
    return chain.GetBranch(name.c_str()) != nullptr;
  }

  bool requireBranches(TChain& chain, const std::vector<std::string>& names)
  {
    bool allPresent = true;
    for (const auto& name : names)
    {
      if (!branchExists(chain, name))
      {
        std::cerr << "ERROR: required branch is missing: " << name << std::endl;
        allPresent = false;
      }
    }
    return allPresent;
  }

  void bindDaughter(TChain& chain,
                    const std::string& prefix,
                    DaughterBranches& daughter)
  {
    chain.SetBranchAddress((prefix + "_track_id").c_str(), &daughter.track_id);
    chain.SetBranchAddress((prefix + "_charge").c_str(), &daughter.charge);
    chain.SetBranchAddress((prefix + "_side").c_str(), &daughter.side);
    chain.SetBranchAddress((prefix + "_npoints").c_str(), &daughter.npoints);
    chain.SetBranchAddress((prefix + "_ntpc_clusters").c_str(), &daughter.ntpc_clusters);
    chain.SetBranchAddress((prefix + "_dedx").c_str(), &daughter.dedx);
    chain.SetBranchAddress((prefix + "_px").c_str(), &daughter.px);
    chain.SetBranchAddress((prefix + "_py").c_str(), &daughter.py);
    chain.SetBranchAddress((prefix + "_pz").c_str(), &daughter.pz);
    chain.SetBranchAddress((prefix + "_pt").c_str(), &daughter.pt);
    chain.SetBranchAddress((prefix + "_eta").c_str(), &daughter.eta);
    chain.SetBranchAddress((prefix + "_fit_chi2").c_str(), &daughter.fit_chi2);
    chain.SetBranchAddress((prefix + "_fit_ndf").c_str(), &daughter.fit_ndf);
    chain.SetBranchAddress((prefix + "_fit_chi2_ndf").c_str(), &daughter.fit_chi2_ndf);

    chain.SetBranchAddress((prefix + "_cluster_index").c_str(), &daughter.cluster_index);
    chain.SetBranchAddress((prefix + "_cluster_side").c_str(), &daughter.cluster_side);
    chain.SetBranchAddress((prefix + "_layer").c_str(), &daughter.layer);

    chain.SetBranchAddress((prefix + "_cluster_x").c_str(), &daughter.cluster_x);
    chain.SetBranchAddress((prefix + "_cluster_y").c_str(), &daughter.cluster_y);
    chain.SetBranchAddress((prefix + "_cluster_z").c_str(), &daughter.cluster_z);
    chain.SetBranchAddress((prefix + "_cluster_r").c_str(), &daughter.cluster_r);
    chain.SetBranchAddress((prefix + "_cluster_phi").c_str(), &daughter.cluster_phi);

    chain.SetBranchAddress((prefix + "_fit_x").c_str(), &daughter.fit_x);
    chain.SetBranchAddress((prefix + "_fit_y").c_str(), &daughter.fit_y);
    chain.SetBranchAddress((prefix + "_fit_z").c_str(), &daughter.fit_z);
    chain.SetBranchAddress((prefix + "_fit_px").c_str(), &daughter.fit_px);
    chain.SetBranchAddress((prefix + "_fit_py").c_str(), &daughter.fit_py);
    chain.SetBranchAddress((prefix + "_fit_pz").c_str(), &daughter.fit_pz);

    chain.SetBranchAddress((prefix + "_residual_r").c_str(), &daughter.residual_r);
    chain.SetBranchAddress((prefix + "_residual_rphi").c_str(), &daughter.residual_rphi);
    chain.SetBranchAddress((prefix + "_residual_z").c_str(), &daughter.residual_z);

    chain.SetBranchAddress((prefix + "_assigned_sigma_r").c_str(), &daughter.assigned_sigma_r);
    chain.SetBranchAddress((prefix + "_assigned_sigma_rphi").c_str(), &daughter.assigned_sigma_rphi);
    chain.SetBranchAddress((prefix + "_assigned_sigma_z").c_str(), &daughter.assigned_sigma_z);

    chain.SetBranchAddress(
        (prefix + "_recommended_for_reference_fit").c_str(),
        &daughter.recommended_for_reference_fit);

    chain.SetBranchAddress(
        (prefix + "_kalman_measurement_chi2").c_str(),
        &daughter.kalman_measurement_chi2);

    chain.SetBranchAddress(
        (prefix + "_kalman_measurement_used").c_str(),
        &daughter.kalman_measurement_used);
  }

  std::vector<std::string> daughterBranchNames(const std::string& prefix)
  {
    return {
        prefix + "_track_id",
        prefix + "_charge",
        prefix + "_side",
        prefix + "_npoints",
        prefix + "_ntpc_clusters",
        prefix + "_dedx",
        prefix + "_px",
        prefix + "_py",
        prefix + "_pz",
        prefix + "_pt",
        prefix + "_eta",
        prefix + "_fit_chi2",
        prefix + "_fit_ndf",
        prefix + "_fit_chi2_ndf",
        prefix + "_cluster_index",
        prefix + "_cluster_side",
        prefix + "_layer",
        prefix + "_cluster_x",
        prefix + "_cluster_y",
        prefix + "_cluster_z",
        prefix + "_cluster_r",
        prefix + "_cluster_phi",
        prefix + "_fit_x",
        prefix + "_fit_y",
        prefix + "_fit_z",
        prefix + "_fit_px",
        prefix + "_fit_py",
        prefix + "_fit_pz",
        prefix + "_residual_r",
        prefix + "_residual_rphi",
        prefix + "_residual_z",
        prefix + "_assigned_sigma_r",
        prefix + "_assigned_sigma_rphi",
        prefix + "_assigned_sigma_z",
        prefix + "_recommended_for_reference_fit",
        prefix + "_kalman_measurement_chi2",
        prefix + "_kalman_measurement_used"};
  }

  std::vector<std::string> categoriesForDaughter(const DaughterBranches& daughter)
  {
    std::vector<std::string> result{"all"};

    if (daughter.side == 0) result.emplace_back("side0");
    if (daughter.side == 1) result.emplace_back("side1");

    if (daughter.charge > 0) result.emplace_back("qplus");
    if (daughter.charge < 0) result.emplace_back("qminus");

    if (daughter.side == 0 && daughter.charge > 0)
      result.emplace_back("side0_qplus");
    if (daughter.side == 0 && daughter.charge < 0)
      result.emplace_back("side0_qminus");
    if (daughter.side == 1 && daughter.charge > 0)
      result.emplace_back("side1_qplus");
    if (daughter.side == 1 && daughter.charge < 0)
      result.emplace_back("side1_qminus");

    return result;
  }

  PairHistSet bookPairHistograms(TDirectory* directory,
                                 const std::string& sampleLabel)
  {
    directory->cd();

    const TString suffix =
        TString::Format(" [%s]", sampleLabel.c_str());

    PairHistSet h;

    h.h_mass = new TH1D(
        "h_mass_Kshort",
        "K^{0}_{S} invariant mass" + suffix +
            ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs",
        400, 0.40, 0.60);

    h.h_v0_pt = new TH1D(
        "h_v0_pt",
        "K^{0}_{S} candidate p_{T}" + suffix +
            ";p_{T}^{K^{0}_{S}} [GeV/c];pairs",
        100, 0.0, 10.0);

    h.h_pca_z = new TH1D(
        "h_pca_z",
        "Secondary-vertex z" + suffix +
            ";PCA z [cm];pairs",
        200, -20.0, 20.0);

    h.h_abs_delta_pca_z = new TH1D(
        "h_abs_delta_pca_z",
        "Daughter PCA-z difference" + suffix +
            ";|PCA_{z,1}-PCA_{z,2}| [cm];pairs",
        200, 0.0, 2.0);

    h.h_pair_dca = new TH1D(
        "h_abs_pairDCA",
        "Daughter-pair DCA" + suffix +
            ";|pair DCA| [cm];pairs",
        200, 0.0, 2.0);

    h.h_decay_radius = new TH1D(
        "h_decay_radius",
        "Transverse decay radius from beam position" + suffix +
            ";R_{decay} [cm];pairs",
        160, 0.0, 80.0);

    h.h_dira = new TH1D(
        "h_dira",
        "Direction angle" + suffix +
            ";DIRA;pairs",
        200, 0.0, 1.0);

    h.h_alpha = new TH1D(
        "h_alpha",
        "Armenteros #alpha" + suffix +
            ";#alpha;pairs",
        240, -1.2, 1.2);

    h.h_qt = new TH1D(
        "h_qT",
        "Armenteros q_{T}" + suffix +
            ";q_{T} [GeV/c];pairs",
        200, 0.0, 0.4);

    h.h_mass_vs_v0_pt = new TH2D(
        "h_mass_Kshort_vs_v0_pt",
        "K^{0}_{S} mass vs candidate p_{T}" + suffix +
            ";p_{T}^{K^{0}_{S}} [GeV/c];m_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
        100, 0.0, 10.0,
        200, 0.40, 0.60);

    h.h_pair_dca_vs_delta_pca_z = new TH2D(
        "h_abs_pairDCA_vs_abs_delta_pca_z",
        "Pair DCA vs daughter PCA-z difference" + suffix +
            ";|PCA_{z,1}-PCA_{z,2}| [cm];|pair DCA| [cm]",
        200, 0.0, 2.0,
        200, 0.0, 2.0);

    h.h_daughter_pt_correlation = new TH2D(
        "h_daughter_pt1_vs_pt2",
        "Daughter p_{T} correlation" + suffix +
            ";p_{T,1} [GeV/c];p_{T,2} [GeV/c]",
        100, 0.0, 5.0,
        100, 0.0, 5.0);

    h.h_dca_xy_correlation = new TH2D(
        "h_dca_xy1_vs_dca_xy2",
        "Daughter DCA_{xy} correlation" + suffix +
            ";DCA_{xy,1} [cm];DCA_{xy,2} [cm]",
        200, -10.0, 10.0,
        200, -10.0, 10.0);

    h.h_dca_z_correlation = new TH2D(
        "h_dca_z1_vs_dca_z2",
        "Daughter DCA_{z} correlation" + suffix +
            ";DCA_{z,1} [cm];DCA_{z,2} [cm]",
        200, -20.0, 20.0,
        200, -20.0, 20.0);

    return h;
  }

  DaughterHistSet bookDaughterHistograms(TDirectory* directory,
                                         const std::string& sampleLabel)
  {
    directory->cd();

    const TString suffix =
        TString::Format(" [%s]", sampleLabel.c_str());

    DaughterHistSet h;

    h.h_pt = new TH1D(
        "h_pt",
        "Pion daughter p_{T}" + suffix +
            ";p_{T} [GeV/c];tracks",
        150, 0.0, 6.0);

    h.h_eta = new TH1D(
        "h_eta",
        "Pion daughter #eta" + suffix +
            ";#eta;tracks",
        120, -2.0, 2.0);

    h.h_phi = new TH1D(
        "h_atan2_px_py",
        "Pion daughter atan2(px,py)" + suffix +
            ";atan2(px,py) [rad];tracks",
        128, -3.2, 3.2);

    h.h_p = new TH1D(
        "h_p",
        "Pion daughter momentum" + suffix +
            ";p [GeV/c];tracks",
        150, 0.0, 10.0);

    h.h_q_over_pt = new TH1D(
        "h_qOverPt",
        "Pion daughter q/p_{T}" + suffix +
            ";q/p_{T} [(GeV/c)^{-1}];tracks",
        200, -5.0, 5.0);

    h.h_npoints = new TH1D(
        "h_npoints",
        "Track point count" + suffix +
            ";n points;tracks",
        64, -0.5, 63.5);

    h.h_ntpc_clusters = new TH1D(
        "h_ntpc_clusters",
        "TPC cluster count" + suffix +
            ";n TPC clusters;tracks",
        64, -0.5, 63.5);

    h.h_fit_chi2_ndf = new TH1D(
        "h_fit_chi2_ndf",
        "Track fit #chi^{2}/ndf" + suffix +
            ";#chi^{2}/ndf;tracks",
        200, 0.0, 50.0);

    h.h_dca_xy_to_primary_vertex = new TH1D(
        "h_dca_xy_to_primary_vertex",
        "Daughter DCA_{xy} to configured primary vertex" + suffix +
            ";DCA_{xy} [cm];tracks",
        240, -12.0, 12.0);

    h.h_dca_z_to_primary_vertex = new TH1D(
        "h_dca_z_to_primary_vertex",
        "Daughter DCA_{z} to configured primary vertex" + suffix +
            ";DCA_{z} [cm];tracks",
        240, -30.0, 30.0);

    h.h_number_of_detail_points = new TH1D(
        "h_number_of_detail_points",
        "Stored daughter-detail point count" + suffix +
            ";vector size;tracks",
        64, -0.5, 63.5);

    h.h_pt_vs_eta = new TH2D(
        "h_pt_vs_eta",
        "Daughter p_{T} vs #eta" + suffix +
            ";#eta;p_{T} [GeV/c]",
        100, -2.0, 2.0,
        120, 0.0, 6.0);

    h.h_dca_xy_vs_phi = new TH2D(
        "h_dca_xy_vs_atan2_px_py",
        "Daughter DCA_{xy} vs atan2(px,py)" + suffix +
            ";atan2(px,py) [rad];DCA_{xy} [cm]",
        128, -3.2, 3.2,
        240, -12.0, 12.0);

    h.h_dca_z_vs_phi = new TH2D(
        "h_dca_z_vs_atan2_px_py",
        "Daughter DCA_{z} vs atan2(px,py)" + suffix +
            ";atan2(px,py) [rad];DCA_{z} [cm]",
        128, -3.2, 3.2,
        240, -30.0, 30.0);

    h.h_dedx_vs_signed_p = new TH2D(
        "h_dedx_vs_signedP",
        "Pion daughter dE/dx vs charge#timesp" + suffix +
            ";qp [GeV/c];dE/dx",
        400, -4.0, 4.0,
        600, 0.0, 6000.0);

    h.h_npoints_vs_ntpc_clusters = new TH2D(
        "h_npoints_vs_ntpc_clusters",
        "Track points vs TPC clusters" + suffix +
            ";n TPC clusters;n points",
        64, -0.5, 63.5,
        64, -0.5, 63.5);

    h.h_fit_chi2_ndf_vs_pt = new TH2D(
        "h_fit_chi2_ndf_vs_pt",
        "Track fit #chi^{2}/ndf vs p_{T}" + suffix +
            ";p_{T} [GeV/c];#chi^{2}/ndf",
        120, 0.0, 6.0,
        200, 0.0, 50.0);

    h.h_primary_dca_xy_vs_phi_vs_pt = new TH3D(
        "h_primary_dca_xy_vs_phi_vs_pt",
        "DCA_{xy} to primary vertex vs #phi and p_{T}" + suffix +
            ";atan2(px,py) [rad];DCA_{xy}^{primary} [cm];p_{T} [GeV/c]",
        96, -3.2, 3.2, 160, -8.0, 8.0, 60, 0.0, 6.0);

    h.h_primary_dca_z_vs_eta_vs_pt = new TH3D(
        "h_primary_dca_z_vs_eta_vs_pt",
        "DCA_{z} to primary vertex vs #eta and p_{T}" + suffix +
            ";#eta;DCA_{z}^{primary} [cm];p_{T} [GeV/c]",
        80, -2.0, 2.0, 180, -18.0, 18.0, 60, 0.0, 6.0);

    h.h_secondary_dca_xy_vs_q_over_pt = new TH2D(
        "h_secondary_dca_xy_vs_q_over_pt",
        "Signed secondary d_{xy} vs q/p_{T}" + suffix +
            ";q/p_{T} [(GeV/c)^{-1}];signed d_{xy}^{secondary} [cm]",
        120, -6.0, 6.0, 200, -1.5, 1.5);

    h.h_secondary_dca_z_vs_q_over_pt = new TH2D(
        "h_secondary_dca_z_vs_q_over_pt",
        "Secondary d_{z} vs q/p_{T}" + suffix +
            ";q/p_{T} [(GeV/c)^{-1}];d_{z}^{secondary} [cm]",
        120, -6.0, 6.0, 160, -0.8, 0.8);

    h.h_secondary_dca_xy_vs_phi_vs_pt = new TH3D(
        "h_secondary_dca_xy_vs_phi_vs_pt",
        "Signed secondary d_{xy} vs #phi and p_{T}" + suffix +
            ";atan2(px,py) [rad];signed d_{xy}^{secondary} [cm];p_{T} [GeV/c]",
        96, -3.2, 3.2, 200, -1.5, 1.5, 60, 0.0, 6.0);

    h.h_secondary_dca_xy_vs_eta_vs_pt = new TH3D(
        "h_secondary_dca_xy_vs_eta_vs_pt",
        "Signed secondary d_{xy} vs #eta and p_{T}" + suffix +
            ";#eta;signed d_{xy}^{secondary} [cm];p_{T} [GeV/c]",
        80, -2.0, 2.0, 200, -1.5, 1.5, 60, 0.0, 6.0);

    h.h_secondary_dca_z_vs_phi_vs_pt = new TH3D(
        "h_secondary_dca_z_vs_phi_vs_pt",
        "Secondary d_{z} vs #phi and p_{T}" + suffix +
            ";atan2(px,py) [rad];d_{z}^{secondary} [cm];p_{T} [GeV/c]",
        96, -3.2, 3.2, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_secondary_dca_z_vs_eta_vs_pt = new TH3D(
        "h_secondary_dca_z_vs_eta_vs_pt",
        "Secondary d_{z} vs #eta and p_{T}" + suffix +
            ";#eta;d_{z}^{secondary} [cm];p_{T} [GeV/c]",
        80, -2.0, 2.0, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.p_secondary_dca_xy_phi_eta_pt0p2 = new TProfile2D(
        "p_secondary_dca_xy_phi_eta_pt0p2",
        "Mean signed secondary d_{xy}, p_{T}>0.2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<d_{xy}^{secondary}> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_secondary_dca_xy_phi_eta_pt1 = new TProfile2D(
        "p_secondary_dca_xy_phi_eta_pt1",
        "Mean signed secondary d_{xy}, p_{T}>1 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<d_{xy}^{secondary}> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_secondary_dca_xy_phi_eta_pt2 = new TProfile2D(
        "p_secondary_dca_xy_phi_eta_pt2",
        "Mean signed secondary d_{xy}, p_{T}>2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<d_{xy}^{secondary}> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");

    h.p_secondary_dca_z_phi_eta_pt0p2 = new TProfile2D(
        "p_secondary_dca_z_phi_eta_pt0p2",
        "Mean secondary d_{z}, p_{T}>0.2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<d_{z}^{secondary}> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_secondary_dca_z_phi_eta_pt1 = new TProfile2D(
        "p_secondary_dca_z_phi_eta_pt1",
        "Mean secondary d_{z}, p_{T}>1 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<d_{z}^{secondary}> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_secondary_dca_z_phi_eta_pt2 = new TProfile2D(
        "p_secondary_dca_z_phi_eta_pt2",
        "Mean secondary d_{z}, p_{T}>2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<d_{z}^{secondary}> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");

    h.h_quality_vs_phi_vs_pt = new TH3D(
        "h_quality_vs_phi_vs_pt",
        "Track quality vs #phi and p_{T}" + suffix +
            ";atan2(px,py) [rad];#chi^{2}/ndf;p_{T} [GeV/c]",
        96, -3.2, 3.2, 160, 0.0, 40.0, 60, 0.0, 6.0);

    h.h_quality_vs_eta_vs_pt = new TH3D(
        "h_quality_vs_eta_vs_pt",
        "Track quality vs #eta and p_{T}" + suffix +
            ";#eta;#chi^{2}/ndf;p_{T} [GeV/c]",
        80, -2.0, 2.0, 160, 0.0, 40.0, 60, 0.0, 6.0);

    h.h_phi_vs_eta_vs_pt = new TH3D(
        "h_phi_vs_eta_vs_pt",
        "Selected-pion phase space" + suffix +
            ";atan2(px,py) [rad];#eta;p_{T} [GeV/c]",
        96, -3.2, 3.2, 80, -2.0, 2.0, 60, 0.0, 6.0);

    h.h_residual_r = new TH1D(
        "h_residual_r",
        "Cluster radial residual" + suffix +
            ";#Deltar [cm];clusters",
        240, -1.2, 1.2);

    h.h_residual_rphi = new TH1D(
        "h_residual_rphi",
        "Cluster r#phi residual" + suffix +
            ";#Deltar#phi [cm];clusters",
        240, -1.2, 1.2);

    h.h_residual_z = new TH1D(
        "h_residual_z",
        "Cluster z residual" + suffix +
            ";#Deltaz [cm];clusters",
        240, -2.0, 2.0);

    h.h_pull_r = new TH1D(
        "h_pull_r",
        "Radial residual pull" + suffix +
            ";#Deltar/#sigma_{r};clusters",
        240, -12.0, 12.0);

    h.h_pull_rphi = new TH1D(
        "h_pull_rphi",
        "r#phi residual pull" + suffix +
            ";#Deltar#phi/#sigma_{r#phi};clusters",
        240, -12.0, 12.0);

    h.h_pull_z = new TH1D(
        "h_pull_z",
        "z residual pull" + suffix +
            ";#Deltaz/#sigma_{z};clusters",
        240, -12.0, 12.0);

    h.h_residual_r_vs_layer = new TH2D(
        "h_residual_r_vs_layer",
        "Radial residual vs layer" + suffix +
            ";TPC layer;#Deltar [cm]",
        48, 6.5, 54.5,
        240, -1.2, 1.2);

    h.h_residual_rphi_vs_layer = new TH2D(
        "h_residual_rphi_vs_layer",
        "r#phi residual vs layer" + suffix +
            ";TPC layer;#Deltar#phi [cm]",
        48, 6.5, 54.5,
        240, -1.2, 1.2);

    h.h_residual_z_vs_layer = new TH2D(
        "h_residual_z_vs_layer",
        "z residual vs layer" + suffix +
            ";TPC layer;#Deltaz [cm]",
        48, 6.5, 54.5,
        240, -2.0, 2.0);

    h.h_pull_r_vs_layer = new TH2D(
        "h_pull_r_vs_layer",
        "Radial residual pull vs layer" + suffix +
            ";TPC layer;#Deltar/#sigma_{r}",
        48, 6.5, 54.5,
        200, -10.0, 10.0);

    h.h_pull_rphi_vs_layer = new TH2D(
        "h_pull_rphi_vs_layer",
        "r#phi residual pull vs layer" + suffix +
            ";TPC layer;#Deltar#phi/#sigma_{r#phi}",
        48, 6.5, 54.5,
        200, -10.0, 10.0);

    h.h_pull_z_vs_layer = new TH2D(
        "h_pull_z_vs_layer",
        "z residual pull vs layer" + suffix +
            ";TPC layer;#Deltaz/#sigma_{z}",
        48, 6.5, 54.5,
        200, -10.0, 10.0);

    h.h_residual_r_vs_cluster_r = new TH2D(
        "h_residual_r_vs_cluster_r",
        "Radial residual vs cluster radius" + suffix +
            ";cluster r [cm];#Deltar [cm]",
        55, 20.0, 75.0, 200, -0.8, 0.8);

    h.h_residual_r_vs_cluster_r_vs_pt = new TH3D(
        "h_residual_r_vs_cluster_r_vs_pt",
        "Radial residual vs cluster radius and p_{T}" + suffix +
            ";cluster r [cm];#Deltar [cm];p_{T} [GeV/c]",
        55, 20.0, 75.0, 200, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_rphi_vs_cluster_r = new TH2D(
        "h_residual_rphi_vs_cluster_r",
        "r#phi residual vs cluster radius" + suffix +
            ";cluster r [cm];#Deltar#phi [cm]",
        100, 25.0, 80.0,
        240, -1.2, 1.2);

    h.h_residual_z_vs_cluster_r = new TH2D(
        "h_residual_z_vs_cluster_r",
        "z residual vs cluster radius" + suffix +
            ";cluster r [cm];#Deltaz [cm]",
        100, 25.0, 80.0,
        240, -2.0, 2.0);

    h.h_residual_rphi_vs_cluster_phi = new TH2D(
        "h_residual_rphi_vs_cluster_phi",
        "r#phi residual vs cluster #phi" + suffix +
            ";cluster #phi [rad];#Deltar#phi [cm]",
        128, -3.2, 3.2,
        240, -1.2, 1.2);

    h.h_residual_z_vs_cluster_phi = new TH2D(
        "h_residual_z_vs_cluster_phi",
        "z residual vs cluster #phi" + suffix +
            ";cluster #phi [rad];#Deltaz [cm]",
        128, -3.2, 3.2,
        240, -2.0, 2.0);

    h.h_residual_z_vs_cluster_z = new TH2D(
        "h_residual_z_vs_cluster_z",
        "z residual vs cluster z" + suffix +
            ";cluster z [cm];#Deltaz [cm]",
        240, -120.0, 120.0,
        240, -2.0, 2.0);
    h.h_residual_rphi_vs_cluster_r_vs_pt = new TH3D(
        "h_residual_rphi_vs_cluster_r_vs_pt",
        "r#phi residual vs cluster radius and p_{T}" + suffix +
            ";cluster r [cm];#Deltar#phi [cm];p_{T} [GeV/c]",
        55, 20.0, 75.0, 200, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_z_vs_cluster_r_vs_pt = new TH3D(
        "h_residual_z_vs_cluster_r_vs_pt",
        "z residual vs cluster radius and p_{T}" + suffix +
            ";cluster r [cm];#Deltaz [cm];p_{T} [GeV/c]",
        55, 20.0, 75.0, 240, -2.0, 2.0, 60, 0.0, 6.0);

    h.h_residual_z_vs_cluster_z_vs_pt = new TH3D(
        "h_residual_z_vs_cluster_z_vs_pt",
        "z residual vs cluster z and p_{T}" + suffix +
            ";cluster z [cm];#Deltaz [cm];p_{T} [GeV/c]",
        240, -120.0, 120.0, 240, -2.0, 2.0, 60, 0.0, 6.0);


    h.h_residual_r_vs_phi_vs_pt = new TH3D(
        "h_residual_r_vs_phi_vs_pt",
        "Radial residual vs cluster #phi and daughter p_{T}" + suffix +
            ";cluster #phi [rad];#Deltar [cm];p_{T} [GeV/c]",
        96, -3.2, 3.2, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_rphi_vs_phi_vs_pt = new TH3D(
        "h_residual_rphi_vs_phi_vs_pt",
        "r#phi residual vs cluster #phi and daughter p_{T}" + suffix +
            ";cluster #phi [rad];#Deltar#phi [cm];p_{T} [GeV/c]",
        96, -3.2, 3.2, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_z_vs_phi_vs_pt = new TH3D(
        "h_residual_z_vs_phi_vs_pt",
        "z residual vs cluster #phi and daughter p_{T}" + suffix +
            ";cluster #phi [rad];#Deltaz [cm];p_{T} [GeV/c]",
        96, -3.2, 3.2, 180, -1.8, 1.8, 60, 0.0, 6.0);

    h.h_residual_r_vs_eta_vs_pt = new TH3D(
        "h_residual_r_vs_eta_vs_pt",
        "Radial residual vs daughter #eta and p_{T}" + suffix +
            ";#eta;#Deltar [cm];p_{T} [GeV/c]",
        80, -2.0, 2.0, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_rphi_vs_eta_vs_pt = new TH3D(
        "h_residual_rphi_vs_eta_vs_pt",
        "r#phi residual vs daughter #eta and p_{T}" + suffix +
            ";#eta;#Deltar#phi [cm];p_{T} [GeV/c]",
        80, -2.0, 2.0, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_z_vs_eta_vs_pt = new TH3D(
        "h_residual_z_vs_eta_vs_pt",
        "z residual vs daughter #eta and p_{T}" + suffix +
            ";#eta;#Deltaz [cm];p_{T} [GeV/c]",
        80, -2.0, 2.0, 180, -1.8, 1.8, 60, 0.0, 6.0);

    h.h_pull_rphi_vs_phi_vs_pt = new TH3D(
        "h_pull_rphi_vs_phi_vs_pt",
        "r#phi pull vs cluster #phi and daughter p_{T}" + suffix +
            ";cluster #phi [rad];#Deltar#phi/#sigma_{r#phi};p_{T} [GeV/c]",
        96, -3.2, 3.2, 160, -8.0, 8.0, 60, 0.0, 6.0);

    h.h_pull_z_vs_eta_vs_pt = new TH3D(
        "h_pull_z_vs_eta_vs_pt",
        "z pull vs daughter #eta and p_{T}" + suffix +
            ";#eta;#Deltaz/#sigma_{z};p_{T} [GeV/c]",
        80, -2.0, 2.0, 160, -8.0, 8.0, 60, 0.0, 6.0);

    h.h_residual_rphi_vs_layer_vs_pt = new TH3D(
        "h_residual_rphi_vs_layer_vs_pt",
        "r#phi residual vs layer and daughter p_{T}" + suffix +
            ";TPC layer;#Deltar#phi [cm];p_{T} [GeV/c]",
        48, 6.5, 54.5, 160, -0.8, 0.8, 60, 0.0, 6.0);

    h.h_residual_z_vs_layer_vs_pt = new TH3D(
        "h_residual_z_vs_layer_vs_pt",
        "z residual vs layer and daughter p_{T}" + suffix +
            ";TPC layer;#Deltaz [cm];p_{T} [GeV/c]",
        48, 6.5, 54.5, 180, -1.8, 1.8, 60, 0.0, 6.0);

    h.h_assigned_sigma_r_vs_layer = new TH2D(
        "h_assigned_sigma_r_vs_layer",
        "Assigned #sigma_{r} vs layer" + suffix +
            ";TPC layer;#sigma_{r} [cm]",
        48, 6.5, 54.5,
        200, 0.0, 1.0);

    h.h_assigned_sigma_rphi_vs_layer = new TH2D(
        "h_assigned_sigma_rphi_vs_layer",
        "Assigned #sigma_{r#phi} vs layer" + suffix +
            ";TPC layer;#sigma_{r#phi} [cm]",
        48, 6.5, 54.5,
        200, 0.0, 0.5);

    h.h_assigned_sigma_z_vs_layer = new TH2D(
        "h_assigned_sigma_z_vs_layer",
        "Assigned #sigma_{z} vs layer" + suffix +
            ";TPC layer;#sigma_{z} [cm]",
        48, 6.5, 54.5,
        200, 0.0, 0.5);

    h.h_cluster_occupancy_layer_phi = new TH2D(
        "h_cluster_occupancy_layer_phi",
        "Selected-pion cluster occupancy" + suffix +
            ";TPC layer;cluster #phi [rad]",
        48, 6.5, 54.5,
        128, -3.2, 3.2);

    h.h_cluster_occupancy_cluster_z_layer = new TH2D(
        "h_cluster_occupancy_cluster_z_layer",
        "Selected-pion cluster occupancy vs z and layer" + suffix +
            ";cluster z [cm];TPC layer",
        240, -120.0, 120.0,
        48, 6.5, 54.5);

    h.p_mean_residual_r_vs_layer = new TProfile(
        "p_mean_residual_r_vs_layer",
        "Mean radial residual vs layer" + suffix +
            ";TPC layer;<#Deltar> [cm]",
        48, 6.5, 54.5, "s");

    h.p_mean_residual_rphi_vs_layer = new TProfile(
        "p_mean_residual_rphi_vs_layer",
        "Mean r#phi residual vs layer" + suffix +
            ";TPC layer;<#Deltar#phi> [cm]",
        48, 6.5, 54.5, "s");

    h.p_mean_residual_z_vs_layer = new TProfile(
        "p_mean_residual_z_vs_layer",
        "Mean z residual vs layer" + suffix +
            ";TPC layer;<#Deltaz> [cm]",
        48, 6.5, 54.5, "s");

    h.p_recommended_fraction_vs_layer = new TProfile(
        "p_recommended_fraction_vs_layer",
        "Reference-fit recommendation fraction vs layer" + suffix +
            ";TPC layer;recommended fraction",
        48, 6.5, 54.5);

    h.p_mean_residual_rphi_layer_phi = new TProfile2D(
        "p_mean_residual_rphi_layer_phi",
        "Mean r#phi residual vs layer and cluster #phi" + suffix +
            ";TPC layer;cluster #phi [rad];<#Deltar#phi> [cm]",
        48, 6.5, 54.5,
        128, -3.2, 3.2,
        "s");

    h.p_mean_residual_z_layer_phi = new TProfile2D(
        "p_mean_residual_z_layer_phi",
        "Mean z residual vs layer and cluster #phi" + suffix +
            ";TPC layer;cluster #phi [rad];<#Deltaz> [cm]",
        48, 6.5, 54.5,
        128, -3.2, 3.2,
        "s");

    h.p_mean_residual_z_cluster_z_layer = new TProfile2D(
        "p_mean_residual_z_cluster_z_layer",
        "Mean z residual vs cluster z and layer" + suffix +
            ";cluster z [cm];TPC layer;<#Deltaz> [cm]",
        240, -120.0, 120.0,
        48, 6.5, 54.5,
        "s");

    h.p_residual_rphi_phi_eta_pt0p2 = new TProfile2D(
        "p_residual_rphi_phi_eta_pt0p2",
        "Mean r#phi residual vs #phi and #eta, p_{T}>0.2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<#Deltar#phi> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_residual_rphi_phi_eta_pt1 = new TProfile2D(
        "p_residual_rphi_phi_eta_pt1",
        "Mean r#phi residual vs #phi and #eta, p_{T}>1 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<#Deltar#phi> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_residual_rphi_phi_eta_pt2 = new TProfile2D(
        "p_residual_rphi_phi_eta_pt2",
        "Mean r#phi residual vs #phi and #eta, p_{T}>2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<#Deltar#phi> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");

    h.p_residual_z_phi_eta_pt0p2 = new TProfile2D(
        "p_residual_z_phi_eta_pt0p2",
        "Mean z residual vs #phi and #eta, p_{T}>0.2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<#Deltaz> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_residual_z_phi_eta_pt1 = new TProfile2D(
        "p_residual_z_phi_eta_pt1",
        "Mean z residual vs #phi and #eta, p_{T}>1 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<#Deltaz> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");
    h.p_residual_z_phi_eta_pt2 = new TProfile2D(
        "p_residual_z_phi_eta_pt2",
        "Mean z residual vs #phi and #eta, p_{T}>2 GeV/c" + suffix +
            ";atan2(px,py) [rad];#eta;<#Deltaz> [cm]",
        96, -3.2, 3.2, 80, -2.0, 2.0, "s");

    h.p_mean_residual_rphi_layer_phi_reference = new TProfile2D(
        "p_mean_residual_rphi_layer_phi_reference",
        "Mean r#phi residual for recommended reference-fit points" + suffix +
            ";TPC layer;cluster #phi [rad];<#Deltar#phi> [cm]",
        48, 6.5, 54.5,
        128, -3.2, 3.2,
        "s");

    h.p_mean_residual_z_layer_phi_reference = new TProfile2D(
        "p_mean_residual_z_layer_phi_reference",
        "Mean z residual for recommended reference-fit points" + suffix +
            ";TPC layer;cluster #phi [rad];<#Deltaz> [cm]",
        48, 6.5, 54.5,
        128, -3.2, 3.2,
        "s");

    h.h_fit_pt_vs_layer = new TH2D(
        "h_fit_pt_vs_layer",
        "Smoothed fitted p_{T} vs layer" + suffix +
            ";TPC layer;fitted p_{T} [GeV/c]",
        48, 6.5, 54.5,
        150, 0.0, 6.0);

    h.h_fit_q_over_pt_vs_layer = new TH2D(
        "h_fit_qOverPt_vs_layer",
        "Smoothed fitted q/p_{T} vs layer" + suffix +
            ";TPC layer;q/p_{T} [(GeV/c)^{-1}]",
        48, 6.5, 54.5,
        200, -5.0, 5.0);

    h.h_fit_pt_over_track_pt_vs_layer = new TH2D(
        "h_fit_pt_over_track_pt_vs_layer",
        "Fitted p_{T} divided by track p_{T}" + suffix +
            ";TPC layer;p_{T}^{fit}/p_{T}^{track}",
        48, 6.5, 54.5,
        200, 0.5, 1.5);

    h.h_fit_pz_vs_layer = new TH2D(
        "h_fit_pz_vs_layer",
        "Smoothed fitted p_{z} vs layer" + suffix +
            ";TPC layer;p_{z}^{fit} [GeV/c]",
        48, 6.5, 54.5,
        200, -5.0, 5.0);

    h.h_kalman_measurement_chi2 = new TH1D(
        "h_kalman_measurement_chi2",
        "Kalman measurement #chi^{2}" + suffix +
            ";measurement #chi^{2};measurements",
        250, 0.0, 250.0);

    h.h_kalman_measurement_chi2_used = new TH1D(
        "h_kalman_measurement_chi2_used",
        "Kalman measurement #chi^{2}, used" + suffix +
            ";measurement #chi^{2};measurements",
        250, 0.0, 250.0);

    h.h_kalman_measurement_chi2_rejected = new TH1D(
        "h_kalman_measurement_chi2_rejected",
        "Kalman measurement #chi^{2}, rejected" + suffix +
            ";measurement #chi^{2};measurements",
        250, 0.0, 250.0);

    h.h_kalman_used_fraction = new TH1D(
        "h_kalman_used_fraction",
        "Fraction of Kalman measurements used" + suffix +
            ";used fraction;tracks",
        101, -0.005, 1.005);

    h.h_kalman_number_of_measurements = new TH1D(
        "h_kalman_number_of_measurements",
        "Number of Kalman measurements" + suffix +
            ";n measurements;tracks",
        64, -0.5, 63.5);

    h.h_kalman_size_minus_detail_points = new TH1D(
        "h_kalman_size_minus_detail_points",
        "Kalman measurement-vector size minus detail-point size" + suffix +
            ";n_{Kalman}-n_{detail};tracks",
        41, -20.5, 20.5);

    h.h_kalman_measurement_chi2_vs_index = new TH2D(
        "h_kalman_measurement_chi2_vs_index",
        "Kalman measurement #chi^{2} vs measurement index" + suffix +
            ";measurement index;measurement #chi^{2}",
        64, -0.5, 63.5,
        250, 0.0, 250.0);

    h.h_kalman_measurement_used_vs_index = new TH2D(
        "h_kalman_measurement_used_vs_index",
        "Kalman measurement-used flag vs measurement index" + suffix +
            ";measurement index;used flag",
        64, -0.5, 63.5,
        2, -0.5, 1.5);

    h.h_kalman_measurement_chi2_vs_pt = new TH2D(
        "h_kalman_measurement_chi2_vs_pt",
        "Kalman measurement #chi^{2} vs daughter p_{T}" + suffix +
            ";p_{T} [GeV/c];measurement #chi^{2}",
        120, 0.0, 6.0,
        250, 0.0, 250.0);

    h.h_kalman_measurement_chi2_vs_index_vs_pt = new TH3D(
        "h_kalman_measurement_chi2_vs_index_vs_pt",
        "Kalman measurement #chi^{2} vs index and daughter p_{T}" + suffix +
            ";measurement index;measurement #chi^{2};p_{T} [GeV/c]",
        64, -0.5, 63.5, 160, 0.0, 80.0, 60, 0.0, 6.0);

    h.h_kalman_measurement_used_vs_index_vs_pt = new TH3D(
        "h_kalman_measurement_used_vs_index_vs_pt",
        "Kalman measurement used flag vs index and daughter p_{T}" + suffix +
            ";measurement index;used flag;p_{T} [GeV/c]",
        64, -0.5, 63.5, 2, -0.5, 1.5, 60, 0.0, 6.0);

    return h;
  }

  std::size_t commonDetailPointCount(const DaughterBranches& d)
  {
    if (!d.layer ||
        !d.cluster_z ||
        !d.cluster_r ||
        !d.cluster_phi ||
        !d.fit_px ||
        !d.fit_py ||
        !d.fit_pz ||
        !d.residual_r ||
        !d.residual_rphi ||
        !d.residual_z ||
        !d.assigned_sigma_r ||
        !d.assigned_sigma_rphi ||
        !d.assigned_sigma_z ||
        !d.recommended_for_reference_fit)
    {
      return 0;
    }

    return std::min({
        d.layer->size(),
        d.cluster_z->size(),
        d.cluster_r->size(),
        d.cluster_phi->size(),
        d.fit_px->size(),
        d.fit_py->size(),
        d.fit_pz->size(),
        d.residual_r->size(),
        d.residual_rphi->size(),
        d.residual_z->size(),
        d.assigned_sigma_r->size(),
        d.assigned_sigma_rphi->size(),
        d.assigned_sigma_z->size(),
        d.recommended_for_reference_fit->size()});
  }

  void fillPairHistograms(PairHistSet& h,
                          const double mass,
                          const double v0Pt,
                          const double pcaZ,
                          const double absDeltaPcaZ,
                          const double absPairDCA,
                          const double decayRadius,
                          const double dira,
                          const double alpha,
                          const double qt,
                          const double pt1,
                          const double pt2,
                          const double dcaXY1,
                          const double dcaXY2,
                          const double dcaZ1,
                          const double dcaZ2)
  {
    if (std::isfinite(mass)) h.h_mass->Fill(mass);
    if (std::isfinite(v0Pt)) h.h_v0_pt->Fill(v0Pt);
    if (std::isfinite(pcaZ)) h.h_pca_z->Fill(pcaZ);
    if (std::isfinite(absDeltaPcaZ))
      h.h_abs_delta_pca_z->Fill(absDeltaPcaZ);
    if (std::isfinite(absPairDCA)) h.h_pair_dca->Fill(absPairDCA);
    if (std::isfinite(decayRadius)) h.h_decay_radius->Fill(decayRadius);
    if (std::isfinite(dira)) h.h_dira->Fill(dira);
    if (std::isfinite(alpha)) h.h_alpha->Fill(alpha);
    if (std::isfinite(qt)) h.h_qt->Fill(qt);

    if (std::isfinite(v0Pt) && std::isfinite(mass))
      h.h_mass_vs_v0_pt->Fill(v0Pt, mass);

    if (std::isfinite(absDeltaPcaZ) && std::isfinite(absPairDCA))
      h.h_pair_dca_vs_delta_pca_z->Fill(absDeltaPcaZ, absPairDCA);

    if (std::isfinite(pt1) && std::isfinite(pt2))
      h.h_daughter_pt_correlation->Fill(pt1, pt2);

    if (std::isfinite(dcaXY1) && std::isfinite(dcaXY2))
      h.h_dca_xy_correlation->Fill(dcaXY1, dcaXY2);

    if (std::isfinite(dcaZ1) && std::isfinite(dcaZ2))
      h.h_dca_z_correlation->Fill(dcaZ1, dcaZ2);
  }

  void fillDaughterHistograms(DaughterHistSet& h,
                              const DaughterBranches& d,
                              const double dcaXY,
                              const double dcaZ,
                              const double secondaryDcaXY,
                              const double secondaryDcaZ,
                              const double maxAbsResidualCm)
  {
    const double pt = d.pt;
    const double p =
        std::sqrt(static_cast<double>(d.px) * d.px +
                  static_cast<double>(d.py) * d.py +
                  static_cast<double>(d.pz) * d.pz);

    // Preserve the convention used in the user's residual macro.
    const double phiPxPy = std::atan2(d.px, d.py);
    const double qOverPt =
        (std::isfinite(pt) && pt > 0.0)
            ? static_cast<double>(d.charge) / pt
            : 0.0;
    const double signedP = static_cast<double>(d.charge) * p;

    if (std::isfinite(pt)) h.h_pt->Fill(pt);
    if (std::isfinite(d.eta)) h.h_eta->Fill(d.eta);
    if (std::isfinite(phiPxPy)) h.h_phi->Fill(phiPxPy);
    if (std::isfinite(p)) h.h_p->Fill(p);
    if (std::isfinite(qOverPt)) h.h_q_over_pt->Fill(qOverPt);
    h.h_npoints->Fill(d.npoints);
    h.h_ntpc_clusters->Fill(d.ntpc_clusters);
    if (std::isfinite(d.fit_chi2_ndf))
      h.h_fit_chi2_ndf->Fill(d.fit_chi2_ndf);
    if (std::isfinite(dcaXY)) h.h_dca_xy_to_primary_vertex->Fill(dcaXY);
    if (std::isfinite(dcaZ)) h.h_dca_z_to_primary_vertex->Fill(dcaZ);

    if (std::isfinite(pt) && std::isfinite(d.eta))
      h.h_pt_vs_eta->Fill(d.eta, pt);

    if (std::isfinite(phiPxPy) && std::isfinite(dcaXY))
      h.h_dca_xy_vs_phi->Fill(phiPxPy, dcaXY);

    if (std::isfinite(phiPxPy) && std::isfinite(dcaZ))
      h.h_dca_z_vs_phi->Fill(phiPxPy, dcaZ);

    if (std::isfinite(signedP) && std::isfinite(d.dedx))
      h.h_dedx_vs_signed_p->Fill(signedP, d.dedx);

    h.h_npoints_vs_ntpc_clusters->Fill(d.ntpc_clusters, d.npoints);

    if (std::isfinite(pt) && std::isfinite(d.fit_chi2_ndf))
      h.h_fit_chi2_ndf_vs_pt->Fill(pt, d.fit_chi2_ndf);

    if (std::isfinite(phiPxPy) && std::isfinite(dcaXY) && std::isfinite(pt))
      h.h_primary_dca_xy_vs_phi_vs_pt->Fill(phiPxPy, dcaXY, pt);

    if (std::isfinite(d.eta) && std::isfinite(dcaZ) && std::isfinite(pt))
      h.h_primary_dca_z_vs_eta_vs_pt->Fill(d.eta, dcaZ, pt);

    if (std::isfinite(qOverPt) && std::isfinite(secondaryDcaXY))
      h.h_secondary_dca_xy_vs_q_over_pt->Fill(qOverPt, secondaryDcaXY);
    if (std::isfinite(qOverPt) && std::isfinite(secondaryDcaZ))
      h.h_secondary_dca_z_vs_q_over_pt->Fill(qOverPt, secondaryDcaZ);

    if (std::isfinite(phiPxPy) && std::isfinite(secondaryDcaXY) && std::isfinite(pt))
      h.h_secondary_dca_xy_vs_phi_vs_pt->Fill(phiPxPy, secondaryDcaXY, pt);
    if (std::isfinite(d.eta) && std::isfinite(secondaryDcaXY) && std::isfinite(pt))
      h.h_secondary_dca_xy_vs_eta_vs_pt->Fill(d.eta, secondaryDcaXY, pt);
    if (std::isfinite(phiPxPy) && std::isfinite(secondaryDcaZ) && std::isfinite(pt))
      h.h_secondary_dca_z_vs_phi_vs_pt->Fill(phiPxPy, secondaryDcaZ, pt);
    if (std::isfinite(d.eta) && std::isfinite(secondaryDcaZ) && std::isfinite(pt))
      h.h_secondary_dca_z_vs_eta_vs_pt->Fill(d.eta, secondaryDcaZ, pt);

    if (std::isfinite(phiPxPy) && std::isfinite(d.eta) &&
        std::isfinite(secondaryDcaXY) && std::isfinite(secondaryDcaZ) &&
        std::isfinite(pt))
    {
      if (pt > 0.2)
      {
        h.p_secondary_dca_xy_phi_eta_pt0p2->Fill(phiPxPy, d.eta, secondaryDcaXY);
        h.p_secondary_dca_z_phi_eta_pt0p2->Fill(phiPxPy, d.eta, secondaryDcaZ);
      }
      if (pt > 1.0)
      {
        h.p_secondary_dca_xy_phi_eta_pt1->Fill(phiPxPy, d.eta, secondaryDcaXY);
        h.p_secondary_dca_z_phi_eta_pt1->Fill(phiPxPy, d.eta, secondaryDcaZ);
      }
      if (pt > 2.0)
      {
        h.p_secondary_dca_xy_phi_eta_pt2->Fill(phiPxPy, d.eta, secondaryDcaXY);
        h.p_secondary_dca_z_phi_eta_pt2->Fill(phiPxPy, d.eta, secondaryDcaZ);
      }
    }

    if (std::isfinite(phiPxPy) && std::isfinite(d.fit_chi2_ndf) && std::isfinite(pt))
      h.h_quality_vs_phi_vs_pt->Fill(phiPxPy, d.fit_chi2_ndf, pt);

    if (std::isfinite(d.eta) && std::isfinite(d.fit_chi2_ndf) && std::isfinite(pt))
      h.h_quality_vs_eta_vs_pt->Fill(d.eta, d.fit_chi2_ndf, pt);

    if (std::isfinite(phiPxPy) && std::isfinite(d.eta) && std::isfinite(pt))
      h.h_phi_vs_eta_vs_pt->Fill(phiPxPy, d.eta, pt);

    const std::size_t nDetail = commonDetailPointCount(d);
    h.h_number_of_detail_points->Fill(static_cast<double>(nDetail));

    for (std::size_t i = 0; i < nDetail; ++i)
    {
      const double layer = static_cast<double>(d.layer->at(i));
      const double clusterR = d.cluster_r->at(i);
      const double clusterPhi = d.cluster_phi->at(i);
      const double clusterZ = d.cluster_z->at(i);

      const double residualR = d.residual_r->at(i);
      const double residualRphi = d.residual_rphi->at(i);
      const double residualZ = d.residual_z->at(i);

      const double sigmaR = d.assigned_sigma_r->at(i);
      const double sigmaRphi = d.assigned_sigma_rphi->at(i);
      const double sigmaZ = d.assigned_sigma_z->at(i);

      const bool recommended =
          d.recommended_for_reference_fit->at(i) != 0U;

      if (std::isfinite(layer) && std::isfinite(clusterPhi))
        h.h_cluster_occupancy_layer_phi->Fill(layer, clusterPhi);

      if (std::isfinite(clusterZ) && std::isfinite(layer))
        h.h_cluster_occupancy_cluster_z_layer->Fill(clusterZ, layer);

      if (std::isfinite(layer))
      {
        h.p_recommended_fraction_vs_layer->Fill(
            layer, recommended ? 1.0 : 0.0);

        if (std::isfinite(sigmaR))
          h.h_assigned_sigma_r_vs_layer->Fill(layer, sigmaR);

        if (std::isfinite(sigmaRphi))
          h.h_assigned_sigma_rphi_vs_layer->Fill(layer, sigmaRphi);

        if (std::isfinite(sigmaZ))
          h.h_assigned_sigma_z_vs_layer->Fill(layer, sigmaZ);
      }

      const bool goodR =
          std::isfinite(residualR) &&
          std::abs(residualR) < maxAbsResidualCm;

      const bool goodRphi =
          std::isfinite(residualRphi) &&
          std::abs(residualRphi) < maxAbsResidualCm;

      const bool goodZ =
          std::isfinite(residualZ) &&
          std::abs(residualZ) < maxAbsResidualCm;

      if (goodR)
      {
        h.h_residual_r->Fill(residualR);
        if (std::isfinite(clusterR))
        {
          h.h_residual_r_vs_cluster_r->Fill(clusterR, residualR);
          if (std::isfinite(pt))
            h.h_residual_r_vs_cluster_r_vs_pt->Fill(clusterR, residualR, pt);
        }

        if (std::isfinite(clusterPhi) && std::isfinite(pt))
          h.h_residual_r_vs_phi_vs_pt->Fill(clusterPhi, residualR, pt);
        if (std::isfinite(d.eta) && std::isfinite(pt))
          h.h_residual_r_vs_eta_vs_pt->Fill(d.eta, residualR, pt);

        if (std::isfinite(layer))
        {
          h.h_residual_r_vs_layer->Fill(layer, residualR);
          h.p_mean_residual_r_vs_layer->Fill(layer, residualR);
        }

        if (std::isfinite(sigmaR) && sigmaR > 0.0)
        {
          const double pull = residualR / sigmaR;
          if (std::isfinite(pull))
          {
            h.h_pull_r->Fill(pull);
            if (std::isfinite(layer))
              h.h_pull_r_vs_layer->Fill(layer, pull);
          }
        }
      }

      if (goodRphi)
      {
        h.h_residual_rphi->Fill(residualRphi);

        if (std::isfinite(clusterPhi) && std::isfinite(pt))
          h.h_residual_rphi_vs_phi_vs_pt->Fill(clusterPhi, residualRphi, pt);
        if (std::isfinite(d.eta) && std::isfinite(pt))
          h.h_residual_rphi_vs_eta_vs_pt->Fill(d.eta, residualRphi, pt);
        if (std::isfinite(layer) && std::isfinite(pt))
          h.h_residual_rphi_vs_layer_vs_pt->Fill(layer, residualRphi, pt);

        if (std::isfinite(layer))
        {
          h.h_residual_rphi_vs_layer->Fill(layer, residualRphi);
          h.p_mean_residual_rphi_vs_layer->Fill(layer, residualRphi);
        }

        if (std::isfinite(clusterR))
        {
          h.h_residual_rphi_vs_cluster_r->Fill(clusterR, residualRphi);
          if (std::isfinite(pt))
            h.h_residual_rphi_vs_cluster_r_vs_pt->Fill(clusterR, residualRphi, pt);
        }

        if (std::isfinite(clusterPhi))
          h.h_residual_rphi_vs_cluster_phi->Fill(
              clusterPhi, residualRphi);

        if (std::isfinite(phiPxPy) && std::isfinite(d.eta) && std::isfinite(pt))
        {
          if (pt > 0.2)
            h.p_residual_rphi_phi_eta_pt0p2->Fill(phiPxPy, d.eta, residualRphi);
          if (pt > 1.0)
            h.p_residual_rphi_phi_eta_pt1->Fill(phiPxPy, d.eta, residualRphi);
          if (pt > 2.0)
            h.p_residual_rphi_phi_eta_pt2->Fill(phiPxPy, d.eta, residualRphi);
        }

        if (std::isfinite(layer) && std::isfinite(clusterPhi))
        {
          h.p_mean_residual_rphi_layer_phi->Fill(
              layer, clusterPhi, residualRphi);

          if (recommended)
          {
            h.p_mean_residual_rphi_layer_phi_reference->Fill(
                layer, clusterPhi, residualRphi);
          }
        }

        if (std::isfinite(sigmaRphi) && sigmaRphi > 0.0)
        {
          const double pull = residualRphi / sigmaRphi;
          if (std::isfinite(pull))
          {
            h.h_pull_rphi->Fill(pull);
            if (std::isfinite(layer))
              h.h_pull_rphi_vs_layer->Fill(layer, pull);
            if (std::isfinite(clusterPhi) && std::isfinite(pt))
              h.h_pull_rphi_vs_phi_vs_pt->Fill(clusterPhi, pull, pt);
          }
        }
      }

      if (goodZ)
      {
        h.h_residual_z->Fill(residualZ);

        if (std::isfinite(clusterPhi) && std::isfinite(pt))
          h.h_residual_z_vs_phi_vs_pt->Fill(clusterPhi, residualZ, pt);
        if (std::isfinite(d.eta) && std::isfinite(pt))
          h.h_residual_z_vs_eta_vs_pt->Fill(d.eta, residualZ, pt);
        if (std::isfinite(layer) && std::isfinite(pt))
          h.h_residual_z_vs_layer_vs_pt->Fill(layer, residualZ, pt);

        if (std::isfinite(layer))
        {
          h.h_residual_z_vs_layer->Fill(layer, residualZ);
          h.p_mean_residual_z_vs_layer->Fill(layer, residualZ);
        }

        if (std::isfinite(clusterR))
        {
          h.h_residual_z_vs_cluster_r->Fill(clusterR, residualZ);
          if (std::isfinite(pt))
            h.h_residual_z_vs_cluster_r_vs_pt->Fill(clusterR, residualZ, pt);
        }

        if (std::isfinite(clusterPhi))
          h.h_residual_z_vs_cluster_phi->Fill(
              clusterPhi, residualZ);

        if (std::isfinite(phiPxPy) && std::isfinite(d.eta) && std::isfinite(pt))
        {
          if (pt > 0.2)
            h.p_residual_z_phi_eta_pt0p2->Fill(phiPxPy, d.eta, residualZ);
          if (pt > 1.0)
            h.p_residual_z_phi_eta_pt1->Fill(phiPxPy, d.eta, residualZ);
          if (pt > 2.0)
            h.p_residual_z_phi_eta_pt2->Fill(phiPxPy, d.eta, residualZ);
        }

        if (std::isfinite(clusterZ))
        {
          h.h_residual_z_vs_cluster_z->Fill(clusterZ, residualZ);
          if (std::isfinite(pt))
            h.h_residual_z_vs_cluster_z_vs_pt->Fill(clusterZ, residualZ, pt);
        }

        if (std::isfinite(layer) && std::isfinite(clusterPhi))
        {
          h.p_mean_residual_z_layer_phi->Fill(
              layer, clusterPhi, residualZ);

          if (recommended)
          {
            h.p_mean_residual_z_layer_phi_reference->Fill(
                layer, clusterPhi, residualZ);
          }
        }

        if (std::isfinite(clusterZ) && std::isfinite(layer))
        {
          h.p_mean_residual_z_cluster_z_layer->Fill(
              clusterZ, layer, residualZ);
        }

        if (std::isfinite(sigmaZ) && sigmaZ > 0.0)
        {
          const double pull = residualZ / sigmaZ;
          if (std::isfinite(pull))
          {
            h.h_pull_z->Fill(pull);
            if (std::isfinite(layer))
              h.h_pull_z_vs_layer->Fill(layer, pull);
            if (std::isfinite(d.eta) && std::isfinite(pt))
              h.h_pull_z_vs_eta_vs_pt->Fill(d.eta, pull, pt);
          }
        }
      }

      const double fitPx = d.fit_px->at(i);
      const double fitPy = d.fit_py->at(i);
      const double fitPz = d.fit_pz->at(i);
      const double fitPt = std::hypot(fitPx, fitPy);

      if (std::isfinite(layer) && std::isfinite(fitPt))
      {
        h.h_fit_pt_vs_layer->Fill(layer, fitPt);

        if (fitPt > 0.0)
        {
          const double fitQOverPt =
              static_cast<double>(d.charge) / fitPt;
          h.h_fit_q_over_pt_vs_layer->Fill(layer, fitQOverPt);
        }

        if (std::isfinite(pt) && pt > 0.0)
          h.h_fit_pt_over_track_pt_vs_layer->Fill(layer, fitPt / pt);
      }

      if (std::isfinite(layer) && std::isfinite(fitPz))
        h.h_fit_pz_vs_layer->Fill(layer, fitPz);
    }

    if (d.kalman_measurement_chi2 && d.kalman_measurement_used)
    {
      const std::size_t nKalman =
          std::min(d.kalman_measurement_chi2->size(),
                   d.kalman_measurement_used->size());

      h.h_kalman_number_of_measurements->Fill(
          static_cast<double>(nKalman));

      h.h_kalman_size_minus_detail_points->Fill(
          static_cast<double>(nKalman) -
          static_cast<double>(nDetail));

      std::size_t nUsed = 0;

      for (std::size_t i = 0; i < nKalman; ++i)
      {
        const double chi2 =
            d.kalman_measurement_chi2->at(i);
        const bool used =
            d.kalman_measurement_used->at(i) != 0U;

        if (used) ++nUsed;

        h.h_kalman_measurement_used_vs_index->Fill(
            static_cast<double>(i), used ? 1.0 : 0.0);

        if (!std::isfinite(chi2))
          continue;

        h.h_kalman_measurement_chi2->Fill(chi2);
        h.h_kalman_measurement_chi2_vs_index->Fill(
            static_cast<double>(i), chi2);

        if (std::isfinite(pt))
        {
          h.h_kalman_measurement_chi2_vs_pt->Fill(pt, chi2);
          h.h_kalman_measurement_chi2_vs_index_vs_pt->Fill(
              static_cast<double>(i), chi2, pt);
          h.h_kalman_measurement_used_vs_index_vs_pt->Fill(
              static_cast<double>(i), used ? 1.0 : 0.0, pt);
        }

        if (used)
          h.h_kalman_measurement_chi2_used->Fill(chi2);
        else
          h.h_kalman_measurement_chi2_rejected->Fill(chi2);
      }

      if (nKalman > 0)
      {
        h.h_kalman_used_fraction->Fill(
            static_cast<double>(nUsed) /
            static_cast<double>(nKalman));
      }
    }
  }

  bool passV0Cut(const V0Cut& cut,
                 const DaughterBranches& daughter1,
                 const DaughterBranches& daughter2,
                 const double pairPt1,
                 const double pairPt2,
                 const int pairNpoints1,
                 const int pairNpoints2,
                 const double pcaX,
                 const double pcaY,
                 const double pcaZ,
                 const double pca1Z,
                 const double pca2Z,
                 const double alpha,
                 const double pairDCA,
                 const double dira,
                 const double quality1,
                 const double quality2,
                 const double beamX,
                 const double beamY)
  {
    const double absDeltaPcaZ = std::abs(pca1Z - pca2Z);
    const double decayRadius =
        std::hypot(pcaX, pcaY);
    const double maxQuality =
        std::max(quality1, quality2);

    return
        std::isfinite(pcaZ) &&
        std::isfinite(absDeltaPcaZ) &&
        std::isfinite(decayRadius) &&
        std::isfinite(alpha) &&
        std::isfinite(pairDCA) &&
        std::isfinite(dira) &&
        std::isfinite(maxQuality) &&
        std::isfinite(pairPt1) &&
        std::isfinite(pairPt2) &&
        std::abs(pcaZ) < cut.maxAbsPcaZ &&
        absDeltaPcaZ < cut.maxAbsDeltaPcaZ &&
        pairPt1 > cut.minDaughterPt &&
        pairPt2 > cut.minDaughterPt &&
        decayRadius > cut.minDecayRadius &&
        std::abs(alpha) < cut.maxAbsAlpha &&
        std::abs(pairDCA) < cut.maxAbsPairDCA &&
        dira > cut.minDIRA &&
        maxQuality < cut.maxQuality &&
        pairNpoints1 > cut.minNpoints &&
        pairNpoints2 > cut.minNpoints;
  }

  bool fillCumulativeCutflow(TH1I& cutflow,
                             const V0Cut& cut,
                             const DaughterBranches& daughter1,
                             const DaughterBranches& daughter2,
                             const double pairPt1,
                             const double pairPt2,
                             const int pairNpoints1,
                             const int pairNpoints2,
                             const bool passPionDedx,
                             const bool passDeltaPhi,
                             const double pcaX,
                             const double pcaY,
                             const double pcaZ,
                             const double pca1Z,
                             const double pca2Z,
                             const double alpha,
                             const double pairDCA,
                             const double dira,
                             const double quality1,
                             const double quality2,
                             const double beamX,
                             const double beamY)
  {
    int bin = 1;
    cutflow.Fill(bin++ - 0.5);

    if (!passPionDedx)
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!passDeltaPhi)
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!(std::isfinite(pairPt1) &&
          std::isfinite(pairPt2) &&
          pairPt1 > cut.minDaughterPt &&
          pairPt2 > cut.minDaughterPt))
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!(std::isfinite(pcaZ) &&
          std::abs(pcaZ) < cut.maxAbsPcaZ))
      return false;
    cutflow.Fill(bin++ - 0.5);

    const double absDeltaPcaZ = std::abs(pca1Z - pca2Z);
    if (!(std::isfinite(absDeltaPcaZ) &&
          absDeltaPcaZ < cut.maxAbsDeltaPcaZ))
      return false;
    cutflow.Fill(bin++ - 0.5);

    const double decayRadius =
        std::hypot(pcaX, pcaY);
    if (!(std::isfinite(decayRadius) &&
          decayRadius > cut.minDecayRadius))
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!(std::isfinite(alpha) &&
          std::abs(alpha) < cut.maxAbsAlpha))
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!(std::isfinite(pairDCA) &&
          std::abs(pairDCA) < cut.maxAbsPairDCA))
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!(std::isfinite(dira) &&
          dira > cut.minDIRA))
      return false;
    cutflow.Fill(bin++ - 0.5);

    const double maxQuality =
        std::max(quality1, quality2);
    if (!(std::isfinite(maxQuality) &&
          maxQuality < cut.maxQuality))
      return false;
    cutflow.Fill(bin++ - 0.5);

    if (!(pairNpoints1 > cut.minNpoints &&
          pairNpoints2 > cut.minNpoints))
      return false;
    cutflow.Fill(bin++ - 0.5);

    return true;
  }
}

void MakeK0sPionResidualKalmanQA_v5(
    const char* inputDir = ".",
    const char* filePattern = "*.root",
    const char* outputDir = "output",
    const char* outputName = "k0s_pion_residual_kalman_qa.root",
    const char* treeName = "pairTree",
    const Long64_t maxEntries = -1,
    const bool requireUnlikeSign = true,
    const double signalMassMin = 0.47,
    const double signalMassMax = 0.53,
    const double detailMassMin = 0.40,
    const double detailMassMax = 0.60,
    const double beamX = 0.158,
    const double beamY = 0.285,
    const double maxAbsResidualCm = 2.0)
{
  TH1::AddDirectory(kTRUE);

  const TString chainPattern =
      TString::Format("%s/%s", inputDir, filePattern);

  TChain chain(treeName);
  const int nFiles = chain.Add(chainPattern);

  if (nFiles <= 0)
  {
    std::cerr << "ERROR: no files matched "
              << chainPattern << std::endl;
    return;
  }

  std::vector<std::string> requiredBranches = {
      "mass_Kshort",
      "candidate_mask",
      "has_kshort_daughter_details",
      "v0_pt",
      "v0_px",
      "v0_py",
      "v0_pz",
      "pca_x",
      "pca_y",
      "pca_z",
      "pca1_x",
      "pca1_y",
      "pca1_z",
      "pca2_x",
      "pca2_y",
      "pca2_z",
      "pairDCA",
      "alpha",
      "qT",
      "cosThetaReco",
      "charge1",
      "charge2",
      "dca_xy1",
      "dca_z1",
      "dca_xy2",
      "dca_z2",
      "quality1",
      "quality2",
      "px1",
      "py1",
      "px2",
      "py2",
      "npoints1",
      "npoints2",
      "dedx_1",
      "dedx_2"};

  const auto daughter1Branches =
      daughterBranchNames("daughter1");
  requiredBranches.insert(
      requiredBranches.end(),
      daughter1Branches.begin(),
      daughter1Branches.end());

  const auto daughter2Branches =
      daughterBranchNames("daughter2");
  requiredBranches.insert(
      requiredBranches.end(),
      daughter2Branches.begin(),
      daughter2Branches.end());

  if (!requireBranches(chain, requiredBranches))
  {
    std::cerr
        << "The input must be produced with "
        << "V0_WRITE_KSHORT_DETAILS=true." << std::endl;
    return;
  }

  Float_t massKshort = 0.F;
  UInt_t candidateMask = 0U;
  Int_t hasKshortDaughterDetails = 0;

  Float_t v0Pt = 0.F;
  Float_t v0Px = 0.F;
  Float_t v0Py = 0.F;
  Float_t v0Pz = 0.F;
  Float_t pcaX = 0.F;
  Float_t pcaY = 0.F;
  Float_t pcaZ = 0.F;
  Float_t pca1X = 0.F;
  Float_t pca1Y = 0.F;
  Float_t pca1Z = 0.F;
  Float_t pca2X = 0.F;
  Float_t pca2Y = 0.F;
  Float_t pca2Z = 0.F;
  Float_t pairDCA = 0.F;
  Float_t alpha = 0.F;
  Float_t qt = 0.F;
  Float_t dira = 0.F;
  Float_t charge1 = 0.F;
  Float_t charge2 = 0.F;
  Float_t dcaXY1 = 0.F;
  Float_t dcaZ1 = 0.F;
  Float_t dcaXY2 = 0.F;
  Float_t dcaZ2 = 0.F;
  Float_t quality1 = 0.F;
  Float_t quality2 = 0.F;

  // Pair-level quantities used by the established V0 analysis cuts.
  Float_t pairPx1 = 0.F;
  Float_t pairPy1 = 0.F;
  Float_t pairPx2 = 0.F;
  Float_t pairPy2 = 0.F;
  Short_t pairNpoints1 = 0;
  Short_t pairNpoints2 = 0;
  Float_t pairDedx1 = 0.F;
  Float_t pairDedx2 = 0.F;

  DaughterBranches daughter1;
  DaughterBranches daughter2;

  chain.SetBranchAddress("mass_Kshort", &massKshort);
  chain.SetBranchAddress("candidate_mask", &candidateMask);
  chain.SetBranchAddress(
      "has_kshort_daughter_details",
      &hasKshortDaughterDetails);

  chain.SetBranchAddress("v0_pt", &v0Pt);
  chain.SetBranchAddress("v0_px", &v0Px);
  chain.SetBranchAddress("v0_py", &v0Py);
  chain.SetBranchAddress("v0_pz", &v0Pz);
  chain.SetBranchAddress("pca_x", &pcaX);
  chain.SetBranchAddress("pca_y", &pcaY);
  chain.SetBranchAddress("pca_z", &pcaZ);
  chain.SetBranchAddress("pca1_x", &pca1X);
  chain.SetBranchAddress("pca1_y", &pca1Y);
  chain.SetBranchAddress("pca1_z", &pca1Z);
  chain.SetBranchAddress("pca2_x", &pca2X);
  chain.SetBranchAddress("pca2_y", &pca2Y);
  chain.SetBranchAddress("pca2_z", &pca2Z);
  chain.SetBranchAddress("pairDCA", &pairDCA);
  chain.SetBranchAddress("alpha", &alpha);
  chain.SetBranchAddress("qT", &qt);
  chain.SetBranchAddress("cosThetaReco", &dira);
  chain.SetBranchAddress("charge1", &charge1);
  chain.SetBranchAddress("charge2", &charge2);
  chain.SetBranchAddress("dca_xy1", &dcaXY1);
  chain.SetBranchAddress("dca_z1", &dcaZ1);
  chain.SetBranchAddress("dca_xy2", &dcaXY2);
  chain.SetBranchAddress("dca_z2", &dcaZ2);
  chain.SetBranchAddress("quality1", &quality1);
  chain.SetBranchAddress("quality2", &quality2);
  chain.SetBranchAddress("px1", &pairPx1);
  chain.SetBranchAddress("py1", &pairPy1);
  chain.SetBranchAddress("px2", &pairPx2);
  chain.SetBranchAddress("py2", &pairPy2);
  chain.SetBranchAddress("npoints1", &pairNpoints1);
  chain.SetBranchAddress("npoints2", &pairNpoints2);
  chain.SetBranchAddress("dedx_1", &pairDedx1);
  chain.SetBranchAddress("dedx_2", &pairDedx2);

  bindDaughter(chain, "daughter1", daughter1);
  bindDaughter(chain, "daughter2", daughter2);

  const std::vector<V0Cut> cuts = {
      {
          "cut03_baseline",
          "|pca_z|<15, pT>0.20, |pca1_z-pca2_z|<0.50, "
          "Rdecay>2, |alpha|<0.99, pairDCA<2.0, "
          "DIRA>0.85, quality<15, npoints>30",
          15.0, 0.50, 0.20, 2.0, 0.99, 2.00, 0.85, 15.0, 30},
      {
          "cut07_pairDCA_5mm",
          "|pca_z|<10, pT>0.20, |pca1_z-pca2_z|<0.20, "
          "Rdecay>2, |alpha|<0.99, pairDCA<0.50, "
          "DIRA>0.95, quality<10, npoints>32",
          10.0, 0.20, 0.20, 2.0, 0.99, 0.50, 0.95, 10.0, 32}};

  const std::vector<std::string> massRegions = {
      "signal_0p47_0p53",
      "outside_0p40_0p47_and_0p53_0p60"};

  const std::vector<std::string> categoryNames = {
      "all",
      "side0",
      "side1",
      "qplus",
      "qminus",
      "side0_qplus",
      "side0_qminus",
      "side1_qplus",
      "side1_qminus"};

  gSystem->mkdir(outputDir, kTRUE);

  const TString outputPath =
      TString::Format("%s/%s", outputDir, outputName);

  std::unique_ptr<TFile> output(
      TFile::Open(outputPath, "RECREATE"));

  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: cannot create output file "
              << outputPath << std::endl;
    return;
  }

  output->cd();

  const TString inputTitle = TString::Format("%s, tree=%s",
                      chainPattern.Data(),
                      treeName);
  TNamed inputInfo("input", inputTitle.Data());
  inputInfo.Write();

  const TString signalTitle = TString::Format("%.6f <= mass_Kshort <= %.6f GeV/c^2",
                      signalMassMin,
                      signalMassMax);
  TNamed signalInfo("signal_mass_region", signalTitle.Data());
  signalInfo.Write();

  const TString sidebandTitle = TString::Format(
          "%.6f <= mass_Kshort < %.6f OR "
          "%.6f < mass_Kshort <= %.6f GeV/c^2",
          detailMassMin,
          signalMassMin,
          signalMassMax,
          detailMassMax);
  TNamed sidebandInfo("outside_mass_region", sidebandTitle.Data());
  sidebandInfo.Write();

  TNamed kalmanIndexWarning(
      "kalman_index_warning",
      "daughter*_kalman_measurement_chi2 and used are stored in Kalman "
      "measurement order. The current pairTree does not store "
      "measurement_original_index, so this macro does not assign those "
      "vectors to detector layer.");
  kalmanIndexWarning.Write();

  TNamed residualDefinition(
      "residual_definition",
      "daughter residuals are cluster minus smoothed fitted state; "
      "they are in-fit residuals, not leave-one-out innovations.");
  residualDefinition.Write();

  TNamed selectionParity(
      "selection_parity",
      "KShort QA applies the same pair-level px/py pT, pair-level npoints, "
      "pion dE/dx < 400, V0 delta-phi requirement, decay radius from the "
      "detector origin, pairDCA, recomputed DIRA, alpha, quality, and PCA-z cuts as "
      "MakeK0sPairHistograms.");
  selectionParity.Write();

  TNamed secondaryVertexDefinition(
      "secondary_vertex_definition",
      "secondary_vertex_x/y/z are the pair midpoint coordinates. In particular, "
      "secondary_vertex_z = 0.5 * (pca1_z + pca2_z).");
  secondaryVertexDefinition.Write();

  TNamed dcaClarification(
      "primary_vs_secondary_dca",
      "h_dca_xy_to_primary_vertex and h_dca_z_to_primary_vertex are daughter ""DCA values to the configured primary vertex. "
      "vertex. They are not constrained by the pair-PCA delta-z cut. "
      "secondary_dca_z is daughter_pca_z minus the secondary-vertex z; therefore "
      "|pca1_z-pca2_z|<0.50 implies |secondary_dca_z|<0.25 cm.");
  dcaClarification.Write();

  TNamed secondaryDcaDefinition(
      "secondary_dca_definition",
      "secondary dca_xy is the signed transverse distance from each daughter "
      "PCA point to the reconstructed KShort midpoint vertex; secondary dca_z "
      "is daughter_pca_z minus pair_vertex_z.");
  secondaryDcaDefinition.Write();

  TParameter<double>("beam_x_cm", beamX).Write();
  TParameter<double>("beam_y_cm", beamY).Write();
  TParameter<double>(
      "maximum_absolute_residual_filled_cm",
      maxAbsResidualCm).Write();

  const double signalWidth =
      signalMassMax - signalMassMin;
  const double sidebandWidth =
      (signalMassMin - detailMassMin) +
      (detailMassMax - signalMassMax);

  TParameter<double>(
      "sideband_scale_to_signal_width",
      sidebandWidth > 0.0
          ? signalWidth / sidebandWidth
          : 0.0)
      .Write();

  std::map<std::string, std::map<std::string, SampleHistSet>>
      histograms;

  std::map<std::string, TH1I*> cutflows;

  for (const auto& cut : cuts)
  {
    TDirectory* cutDirectory =
        output->mkdir(cut.name.c_str());
    cutDirectory->cd();

    TNamed cutInfo("selection", cut.description.c_str());
    cutInfo.Write();

    TH1I* cutflow = new TH1I(
        "h_cutflow",
        ("Cumulative cut flow [" + cut.name + "]"
         ";cut;pair count")
            .c_str(),
        12, 0.0, 12.0);

    const std::vector<std::string> cutLabels = {
        "details+KShort+charge",
        "pion dE/dx",
        "V0 delta phi",
        "daughter pT",
        "|pca z|",
        "|delta pca z|",
        "decay radius",
        "|alpha|",
        "|pair DCA|",
        "DIRA",
        "quality",
        "npoints"};

    for (std::size_t index = 0;
         index < cutLabels.size();
         ++index)
    {
      cutflow->GetXaxis()->SetBinLabel(
          static_cast<int>(index + 1),
          cutLabels[index].c_str());
    }

    cutflows[cut.name] = cutflow;

    for (const auto& region : massRegions)
    {
      TDirectory* regionDirectory =
          cutDirectory->mkdir(region.c_str());

      SampleHistSet sample;

      TDirectory* pairDirectory =
          regionDirectory->mkdir("pair");
      sample.pair = bookPairHistograms(
          pairDirectory,
          cut.name + "/" + region);

      TDirectory* daughtersDirectory =
          regionDirectory->mkdir("daughters");

      for (const auto& category : categoryNames)
      {
        TDirectory* categoryDirectory =
            daughtersDirectory->mkdir(category.c_str());

        sample.daughters.emplace(
            category,
            bookDaughterHistograms(
                categoryDirectory,
                cut.name + "/" + region + "/" + category));
      }

      histograms[cut.name][region] =
          std::move(sample);
    }
  }

  const Long64_t totalEntries =
      chain.GetEntries();

  const Long64_t entriesToProcess =
      maxEntries >= 0
          ? std::min(totalEntries, maxEntries)
          : totalEntries;

  std::cout
      << "Added " << nFiles
      << " files, total pairTree entries = "
      << totalEntries
      << ", processing = "
      << entriesToProcess
      << std::endl;

  Long64_t numberWithDetails = 0;
  Long64_t numberInSignal = 0;
  Long64_t numberInSideband = 0;

  for (Long64_t entry = 0;
       entry < entriesToProcess;
       ++entry)
  {
    chain.GetEntry(entry);

    if (entry % 100000 == 0)
    {
      std::cout
          << "Processing " << entry
          << " / " << entriesToProcess
          << std::endl;
    }

    const bool hasDetails =
        hasKshortDaughterDetails != 0;

    const bool hasKshortMask =
        (candidateMask & kCandidateKShort) != 0U;

    const bool unlikeSign =
        charge1 * charge2 < 0.0F;

    if (!hasDetails || !hasKshortMask)
      continue;

    if (requireUnlikeSign && !unlikeSign)
      continue;

    ++numberWithDetails;

    const bool inSignal =
        std::isfinite(massKshort) &&
        massKshort >= signalMassMin &&
        massKshort <= signalMassMax;

    const bool inSideband =
        std::isfinite(massKshort) &&
        ((massKshort >= detailMassMin &&
          massKshort < signalMassMin) ||
         (massKshort > signalMassMax &&
          massKshort <= detailMassMax));

    if (!inSignal && !inSideband)
      continue;

    if (inSignal) ++numberInSignal;
    if (inSideband) ++numberInSideband;

    // Exact pair-level quantities used in MakeK0sPairHistograms.
    const double pairPt1 = std::hypot(pairPx1, pairPy1);
    const double pairPt2 = std::hypot(pairPx2, pairPy2);

    const double phiPos =
        std::atan2(charge1 > 0.F ? pairPy1 : pairPy2,
                   charge1 > 0.F ? pairPx1 : pairPx2);
    const double phiNeg =
        std::atan2(charge1 > 0.F ? pairPy2 : pairPy1,
                   charge1 > 0.F ? pairPx2 : pairPx1);

    double deltaPhi = phiPos - phiNeg;
    while (deltaPhi > TMath::Pi()) deltaPhi -= 2.0 * TMath::Pi();
    while (deltaPhi <= -TMath::Pi()) deltaPhi += 2.0 * TMath::Pi();

    const bool passV0DeltaPhi =
        deltaPhi >= 0.8 - 0.4 * (v0Pt < 2.0 ? v0Pt : 2.0);

    // Same pion PID used in the established KShort histogram macro.
    const bool passPionDedx =
        std::isfinite(pairDedx1) &&
        std::isfinite(pairDedx2) &&
        pairDedx1 < 400.0 &&
        pairDedx2 < 400.0;

    const double v0Momentum =
        std::sqrt(v0Px * v0Px + v0Py * v0Py + v0Pz * v0Pz);
    const double flightLength =
        std::sqrt(pcaX * pcaX + pcaY * pcaY + pcaZ * pcaZ);

    const double exactDira =
        (v0Momentum > 0.0 && flightLength > 0.0)
            ? (v0Px * pcaX + v0Py * pcaY + v0Pz * pcaZ) /
                  (v0Momentum * flightLength)
            : -2.0;

    const double absDeltaPcaZ =
        std::abs(pca1Z - pca2Z);

    const double absPairDCA =
        std::abs(pairDCA);

    // Match MakeK0sPairHistograms exactly.
    const double decayRadius =
        std::hypot(pcaX, pcaY);

    // Daughter offset from the reconstructed secondary KShort vertex.
    // The pair vertex is the midpoint of the two daughter PCA positions.
    const auto signedSecondaryDcaXY =
        [](const DaughterBranches& daughter,
           const double daughterPcaX,
           const double daughterPcaY,
           const double vertexX,
           const double vertexY)
    {
      const double dx = daughterPcaX - vertexX;
      const double dy = daughterPcaY - vertexY;
      const double distance = std::hypot(dx, dy);
      const double cross =
          static_cast<double>(daughter.px) * dy -
          static_cast<double>(daughter.py) * dx;
      return cross >= 0.0 ? distance : -distance;
    };

    const double secondaryDcaXY1 =
        signedSecondaryDcaXY(daughter1, pca1X, pca1Y, pcaX, pcaY);
    const double secondaryDcaXY2 =
        signedSecondaryDcaXY(daughter2, pca2X, pca2Y, pcaX, pcaY);
    const double secondaryDcaZ1 = pca1Z - pcaZ;
    const double secondaryDcaZ2 = pca2Z - pcaZ;

    for (const auto& cut : cuts)
    {
      TH1I* cutflow =
          cutflows.at(cut.name);

      const bool passesCutflow =
          fillCumulativeCutflow(
              *cutflow,
              cut,
              daughter1,
              daughter2,
              pairPt1,
              pairPt2,
              static_cast<int>(pairNpoints1),
              static_cast<int>(pairNpoints2),
              passPionDedx,
              passV0DeltaPhi,
              pcaX,
              pcaY,
              pcaZ,
              pca1Z,
              pca2Z,
              alpha,
              pairDCA,
              exactDira,
              quality1,
              quality2,
              beamX,
              beamY);

      if (!passesCutflow)
        continue;

      // These two requirements were present in the established V0 mass
      // analysis and were missing in the previous pion-QA macro.
      if (!passPionDedx || !passV0DeltaPhi)
        continue;

      // Defensive consistency check against the compact selection function.
      if (!passV0Cut(
              cut,
              daughter1,
              daughter2,
              pairPt1,
              pairPt2,
              static_cast<int>(pairNpoints1),
              static_cast<int>(pairNpoints2),
              pcaX,
              pcaY,
              pcaZ,
              pca1Z,
              pca2Z,
              alpha,
              pairDCA,
              exactDira,
              quality1,
              quality2,
              beamX,
              beamY))
      {
        continue;
      }

      const std::string region =
          inSignal
              ? "signal_0p47_0p53"
              : "outside_0p40_0p47_and_0p53_0p60";

      SampleHistSet& sample =
          histograms.at(cut.name).at(region);

      fillPairHistograms(
          sample.pair,
          massKshort,
          v0Pt,
          pcaZ,
          absDeltaPcaZ,
          absPairDCA,
          decayRadius,
          exactDira,
          alpha,
          qt,
          pairPt1,
          pairPt2,
          dcaXY1,
          dcaXY2,
          dcaZ1,
          dcaZ2);

      for (const auto& category :
           categoriesForDaughter(daughter1))
      {
        fillDaughterHistograms(
            sample.daughters.at(category),
            daughter1,
            dcaXY1,
            dcaZ1,
            secondaryDcaXY1,
            secondaryDcaZ1,
            maxAbsResidualCm);
      }

      for (const auto& category :
           categoriesForDaughter(daughter2))
      {
        fillDaughterHistograms(
            sample.daughters.at(category),
            daughter2,
            dcaXY2,
            dcaZ2,
            secondaryDcaXY2,
            secondaryDcaZ2,
            maxAbsResidualCm);
      }
    }
  }

  output->cd();

  TParameter<Long64_t>(
      "pairs_with_kshort_daughter_details",
      numberWithDetails)
      .Write();

  TParameter<Long64_t>(
      "pairs_in_signal_mass_region_before_cut03_cut07",
      numberInSignal)
      .Write();

  TParameter<Long64_t>(
      "pairs_in_outside_mass_region_before_cut03_cut07",
      numberInSideband)
      .Write();

  output->Write();
  output->Close();

  std::cout
      << "Wrote K0S pion residual/Kalman QA to "
      << outputPath
      << std::endl;
}
