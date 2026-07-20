// MakeResidualHistograms.C
//
// Run, for example:
//   root -l -b -q 'MakeResidualHistograms.C("/path/to/input","*.root","output","qa.root")'
//
// Arguments:
//   inputDir      directory containing ROOT files
//   filePattern   wildcard accepted by TChain::Add, e.g. "*.root" or "residuals_*.root"
//   outputDir     output directory; it is created if needed
//   outputName    output ROOT filename
//   treeName      tree name, default "residuals"
//
// The TChain is traversed exactly once. Cluster-vector histograms are filled
// inside the same entry loop.

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  struct HistSet
  {
    // Track-level plots.
    TH2D* h_side_pca_z = nullptr;
    TH3D* h_side_pca_z_pt = nullptr;

    TH2D* h_rDCAzero_phi = nullptr;
    TH3D* h_rDCAzero_phi_pt = nullptr;

    TH2D* h_zDCA_phi = nullptr;
    TH3D* h_zDCA_phi_pt = nullptr;

    TH2D* h_zDCA_eta = nullptr;
    TH3D* h_zDCA_eta_pt = nullptr;

    TH2D* h_rDCA_qOverPt = nullptr;
    TH3D* h_rDCA_qOverPt_pt = nullptr;

    // Distortion phase-space QA with |rDCA_zero| < 3 cm, |pca_z| < 10 cm,
    // and ntpc_clusters > 30. Existing side/charge categories allow direct
    // positive-to-negative comparisons for each TPC side.
    TH2D* h_pt_phi_rdca3 = nullptr;
    TH2D* h_pt_eta_rdca3 = nullptr;
    TH2D* h_rDCAzero_phi_rdca3 = nullptr;
    TH2D* h_rDCAzero_eta_rdca3 = nullptr;
    TH3D* h_rDCAzero_qOverPt_phi_rdca3 = nullptr;
    TH3D* h_rDCAzero_qOverPt_eta_rdca3 = nullptr;
    TProfile2D* p_mean_rDCAzero_phi_eta_rdca3 = nullptr;
    TProfile2D* p_mean_qOverPt_phi_eta_rdca3 = nullptr;

    TH2D* h_dedx_signedP = nullptr;
    TH3D* h_dedx_signedP_pt = nullptr;

    // Cluster-level plots.
    TH2D* h_residual_rphi_cluster_r = nullptr;
    TH3D* h_residual_rphi_cluster_r_pt = nullptr;

    TH2D* h_residual_z_cluster_r = nullptr;
    TH3D* h_residual_z_cluster_r_pt = nullptr;

    TH2D* h_residual_rphi_layer = nullptr;
    TH3D* h_residual_rphi_layer_pt = nullptr;

    TH2D* h_residual_z_layer = nullptr;
    TH3D* h_residual_z_layer_pt = nullptr;

    TH3D* h_pt_phi_eta_rdca3 = nullptr;

    // Mean residual maps versus detector layer and cluster phi.
    // Two pT selections are stored separately: pT > 1 and pT > 2 GeV/c.
    // Option "s" makes profile-bin errors equal to the residual spread
    // (standard deviation), rather than the error on the mean.
    TProfile2D* h_mean_residual_rphi_layer_phi_pt1 = nullptr;
    TProfile2D* h_mean_residual_rphi_layer_phi_pt2 = nullptr;
    TProfile2D* h_mean_residual_z_layer_phi_pt1 = nullptr;
    TProfile2D* h_mean_residual_z_layer_phi_pt2 = nullptr;

    // Additional loose-selection maps. These do not use the baseline
    // ntpc_clusters >= 30 requirement applied to the existing plots.
    TProfile2D* h_mean_residual_rphi_layer_phi_pt0 = nullptr;
    TProfile2D* h_mean_residual_z_layer_phi_pt0 = nullptr;

    // Mean z-residual maps versus cluster z and detector layer.
    TProfile2D* h_mean_residual_z_cluster_z_layer_pt1 = nullptr;
    TProfile2D* h_mean_residual_z_cluster_z_layer_pt2 = nullptr;
  };

  HistSet bookHistograms(TDirectory* dir, const std::string& category)
  {
    dir->cd();
    HistSet h;

    const TString suffix = TString::Format(" [%s]", category.c_str());

    h.h_side_pca_z =
      new TH2D("h_side_vs_pca_z",
               "TPC side vs PCA z" + suffix + ";pca_z [cm];side",
               200, -100., 100., 2, 0., 2.);
    h.h_side_pca_z_pt =
      new TH3D("h_side_vs_pca_z_vs_pt",
               "TPC side vs PCA z vs p_{T}" + suffix + ";pca_z [cm];side;p_{T} [GeV/c]",
               200, -100., 100., 2, 0., 2., 100, 0., 10.);

    // The user examples use atan2(px,py), so this macro preserves that convention.
    h.h_rDCAzero_phi =
      new TH2D("h_rDCAzero_vs_atan2_px_py",
               "rDCA_{zero} vs atan2(px,py)" + suffix + ";atan2(px,py);rDCA_{zero} [cm]",
               128, -TMath::Pi(), TMath::Pi(), 200, -10., 10.);
    h.h_rDCAzero_phi_pt =
      new TH3D("h_rDCAzero_vs_atan2_px_py_vs_pt",
               "rDCA_{zero} vs atan2(px,py) vs p_{T}" + suffix +
               ";atan2(px,py);rDCA_{zero} [cm];p_{T} [GeV/c]",
               128, -TMath::Pi(), TMath::Pi(), 200, -10., 10., 100, 0., 10.);

    h.h_zDCA_phi =
      new TH2D("h_zDCA_vs_atan2_px_py",
               "zDCA vs atan2(px,py)" + suffix + ";atan2(px,py);zDCA [cm]",
               128, -TMath::Pi(), TMath::Pi(), 200, -10., 10.);
    h.h_zDCA_phi_pt =
      new TH3D("h_zDCA_vs_atan2_px_py_vs_pt",
               "zDCA vs atan2(px,py) vs p_{T}" + suffix +
               ";atan2(px,py);zDCA [cm];p_{T} [GeV/c]",
               128, -TMath::Pi(), TMath::Pi(), 200, -10., 10., 100, 0., 10.);

    h.h_zDCA_eta =
      new TH2D("h_zDCA_vs_eta",
               "zDCA vs #eta" + suffix + ";#eta;zDCA [cm]",
               100, -2., 2., 200, -10., 10.);
    h.h_zDCA_eta_pt =
      new TH3D("h_zDCA_vs_eta_vs_pt",
               "zDCA vs #eta vs p_{T}" + suffix + ";#eta;zDCA [cm];p_{T} [GeV/c]",
               100, -2., 2., 200, -10., 10., 100, 0., 10.);

    h.h_rDCA_qOverPt =
      new TH2D("h_rDCA_vs_qOverPt",
               "rDCA vs q/p_{T}" + suffix + ";q/p_{T} [(GeV/c)^{-1}];rDCA [cm]",
               100, -5., 5., 200, -10., 10.);
    h.h_rDCA_qOverPt_pt =
      new TH3D("h_rDCA_vs_qOverPt_vs_pt",
               "rDCA vs q/p_{T} vs p_{T}" + suffix +
               ";q/p_{T} [(GeV/c)^{-1}];rDCA [cm];p_{T} [GeV/c]",
               100, -5., 5., 200, -10., 10., 100, 0., 10.);

    h.h_pt_phi_rdca3 =
      new TH2D("h_pt_vs_atan2_px_py_rdca3",
               "p_{T} vs atan2(px,py), |rDCA_{zero}|<3 cm" + suffix +
               ";atan2(px,py);p_{T} [GeV/c]",
               128, -TMath::Pi(), TMath::Pi(), 100, 0., 10.);

    h.h_pt_eta_rdca3 =
      new TH2D("h_pt_vs_eta_rdca3",
               "p_{T} vs #eta, |rDCA_{zero}|<3 cm" + suffix +
               ";#eta;p_{T} [GeV/c]",
               100, -2., 2., 100, 0., 10.);

    h.h_rDCAzero_phi_rdca3 =
      new TH2D("h_rDCAzero_vs_atan2_px_py_rdca3",
               "rDCA_{zero} vs atan2(px,py), |rDCA_{zero}|<3 cm" + suffix +
               ";atan2(px,py);rDCA_{zero} [cm]",
               128, -TMath::Pi(), TMath::Pi(), 120, -3., 3.);

    h.h_rDCAzero_eta_rdca3 =
      new TH2D("h_rDCAzero_vs_eta_rdca3",
               "rDCA_{zero} vs #eta, |rDCA_{zero}|<3 cm" + suffix +
               ";#eta;rDCA_{zero} [cm]",
               100, -2., 2., 120, -3., 3.);

    h.h_rDCAzero_qOverPt_phi_rdca3 =
      new TH3D("h_rDCAzero_vs_qOverPt_vs_atan2_px_py_rdca3",
               "rDCA_{zero} vs q/p_{T} vs atan2(px,py), |rDCA_{zero}|<3 cm" + suffix +
               ";atan2(px,py);q/p_{T} [(GeV/c)^{-1}];rDCA_{zero} [cm]",
               128, -TMath::Pi(), TMath::Pi(), 100, -5., 5., 120, -3., 3.);

    h.h_rDCAzero_qOverPt_eta_rdca3 =
      new TH3D("h_rDCAzero_vs_qOverPt_vs_eta_rdca3",
               "rDCA_{zero} vs q/p_{T} vs #eta, |rDCA_{zero}|<3 cm" + suffix +
               ";#eta;q/p_{T} [(GeV/c)^{-1}];rDCA_{zero} [cm]",
               100, -2., 2., 100, -5., 5., 120, -3., 3.);

    h.p_mean_rDCAzero_phi_eta_rdca3 =
      new TProfile2D("p_mean_rDCAzero_vs_atan2_px_py_eta_rdca3",
                     "Mean rDCA_{zero} vs atan2(px,py) and #eta, |rDCA_{zero}|<3 cm" + suffix +
                     ";atan2(px,py);#eta;<rDCA_{zero}> [cm]",
                     128, -TMath::Pi(), TMath::Pi(),
                     100, -2., 2.,
                     "s");

    h.p_mean_qOverPt_phi_eta_rdca3 =
      new TProfile2D("p_mean_qOverPt_vs_atan2_px_py_eta_rdca3",
                     "Mean q/p_{T} vs atan2(px,py) and #eta, |rDCA_{zero}|<3 cm" + suffix +
                     ";atan2(px,py);#eta;<q/p_{T}> [(GeV/c)^{-1}]",
                     128, -TMath::Pi(), TMath::Pi(),
                     100, -2., 2.,
                     "s");

    h.h_pt_phi_eta_rdca3 =
      new TH3D(
        "h_pt_vs_atan2_px_py_vs_eta_rdca3",
        "p_{T} vs atan2(px,py) and #eta, |rDCA_{zero}|<3"
        + suffix
        + ";atan2(px,py) [rad];#eta;p_{T} [GeV/c]",
        120, -TMath::Pi(), TMath::Pi(),
        50, -2., 2.,
        50, 0., 5.
      );

    h.h_dedx_signedP =
      new TH2D("h_dedx_vs_signedP",
               "dE/dx vs charge#timesp" + suffix +
               ";qp [GeV/c];dE/dx",
               400, -2., 2., 1000, 0., 6000.);
    h.h_dedx_signedP_pt =
      new TH3D("h_dedx_vs_signedP_vs_pt",
               "dE/dx vs charge#timesp vs p_{T}" + suffix +
               ";qp [GeV/c];dE/dx;p_{T} [GeV/c]",
               400, -2., 2., 1000, 0., 6000., 100, 0., 10.);

    h.h_residual_rphi_cluster_r =
      new TH2D("h_residual_rphi_vs_cluster_r",
               "r#phi residual vs cluster radius" + suffix +
               ";cluster r [cm];r#phi residual [cm]",
               55, 20., 75., 200, -0.8, 0.8);
    h.h_residual_rphi_cluster_r_pt =
      new TH3D("h_residual_rphi_vs_cluster_r_vs_pt",
               "r#phi residual vs cluster radius vs p_{T}" + suffix +
               ";cluster r [cm];r#phi residual [cm];p_{T} [GeV/c]",
               55, 20., 75., 200, -0.8, 0.8, 100, 0., 10.);

    h.h_residual_z_cluster_r =
      new TH2D("h_residual_z_vs_cluster_r",
               "z residual vs cluster radius" + suffix +
               ";cluster r [cm];z residual [cm]",
               55, 20., 75., 200, -2., 2.);
    h.h_residual_z_cluster_r_pt =
      new TH3D("h_residual_z_vs_cluster_r_vs_pt",
               "z residual vs cluster radius vs p_{T}" + suffix +
               ";cluster r [cm];z residual [cm];p_{T} [GeV/c]",
               55, 20., 75., 200, -2., 2., 100, 0., 10.);

    h.h_residual_rphi_layer =
      new TH2D("h_residual_rphi_vs_layer",
               "r#phi residual vs layer" + suffix +
               ";TPC layer;r#phi residual [cm]",
               48, 6.5, 54.5, 200, -0.8, 0.8);
    h.h_residual_rphi_layer_pt =
      new TH3D("h_residual_rphi_vs_layer_vs_pt",
               "r#phi residual vs layer vs p_{T}" + suffix +
               ";TPC layer;r#phi residual [cm];p_{T} [GeV/c]",
               48, 6.5, 54.5, 200, -0.8, 0.8, 100, 0., 10.);

    h.h_residual_z_layer =
      new TH2D("h_residual_z_vs_layer",
               "z residual vs layer" + suffix +
               ";TPC layer;z residual [cm]",
               48, 6.5, 54.5, 200, -2., 2.);
    h.h_residual_z_layer_pt =
      new TH3D("h_residual_z_vs_layer_vs_pt",
               "z residual vs layer vs p_{T}" + suffix +
               ";TPC layer;z residual [cm];p_{T} [GeV/c]",
               48, 6.5, 54.5, 200, -2., 2., 100, 0., 10.);

    h.h_mean_residual_rphi_layer_phi_pt1 =
      new TProfile2D("h_mean_residual_rphi_vs_layer_phi_pt1",
                     "Mean r#phi residual vs layer and cluster #phi, p_{T}>1 GeV/c" + suffix +
                     ";TPC layer;cluster #phi [rad];<#Deltar#phi> [cm]",
                     48, 6.5, 54.5,
                     128, -TMath::Pi(), TMath::Pi(),
                     "s");

    h.h_mean_residual_rphi_layer_phi_pt2 =
      new TProfile2D("h_mean_residual_rphi_vs_layer_phi_pt2",
                     "Mean r#phi residual vs layer and cluster #phi, p_{T}>2 GeV/c" + suffix +
                     ";TPC layer;cluster #phi [rad];<#Deltar#phi> [cm]",
                     48, 6.5, 54.5,
                     128, -TMath::Pi(), TMath::Pi(),
                     "s");

    h.h_mean_residual_z_layer_phi_pt1 =
      new TProfile2D("h_mean_residual_z_vs_layer_phi_pt1",
                     "Mean z residual vs layer and cluster #phi, p_{T}>1 GeV/c" + suffix +
                     ";TPC layer;cluster #phi [rad];<#Deltaz> [cm]",
                     48, 6.5, 54.5,
                     128, -TMath::Pi(), TMath::Pi(),
                     "s");

    h.h_mean_residual_z_layer_phi_pt2 =
      new TProfile2D("h_mean_residual_z_vs_layer_phi_pt2",
                     "Mean z residual vs layer and cluster #phi, p_{T}>2 GeV/c" + suffix +
                     ";TPC layer;cluster #phi [rad];<#Deltaz> [cm]",
                     48, 6.5, 54.5,
                     128, -TMath::Pi(), TMath::Pi(),
                     "s");

    h.h_mean_residual_rphi_layer_phi_pt0 =
      new TProfile2D("h_mean_residual_rphi_vs_layer_phi_pt0",
                     "Mean r#phi residual vs layer and cluster #phi, p_{T}>0 GeV/c" + suffix +
                     ";TPC layer;cluster #phi [rad];<#Deltar#phi> [cm]",
                     48, 6.5, 54.5,
                     128, -TMath::Pi(), TMath::Pi(),
                     "s");

    h.h_mean_residual_z_layer_phi_pt0 =
      new TProfile2D("h_mean_residual_z_vs_layer_phi_pt0",
                     "Mean z residual vs layer and cluster #phi, p_{T}>0 GeV/c" + suffix +
                     ";TPC layer;cluster #phi [rad];<#Deltaz> [cm]",
                     48, 6.5, 54.5,
                     128, -TMath::Pi(), TMath::Pi(),
                     "s");

    h.h_mean_residual_z_cluster_z_layer_pt1 =
      new TProfile2D("h_mean_residual_z_vs_cluster_z_layer_pt1",
                     "Mean z residual vs cluster z and layer, p_{T}>1 GeV/c" + suffix +
                     ";cluster z [cm];TPC layer;<#Deltaz> [cm]",
                     240, -120., 120.,
                     48, 6.5, 54.5,
                     "s");

    h.h_mean_residual_z_cluster_z_layer_pt2 =
      new TProfile2D("h_mean_residual_z_vs_cluster_z_layer_pt2",
                     "Mean z residual vs cluster z and layer, p_{T}>2 GeV/c" + suffix +
                     ";cluster z [cm];TPC layer;<#Deltaz> [cm]",
                     240, -120., 120.,
                     48, 6.5, 54.5,
                     "s");

    return h;
  }

  std::vector<std::string> categoriesForTrack(const int side, const double charge)
  {
    std::vector<std::string> categories = {"all"};

    if (side == 0) categories.emplace_back("side0");
    if (side == 1) categories.emplace_back("side1");

    if (charge > 0) categories.emplace_back("qplus");
    if (charge < 0) categories.emplace_back("qminus");

    if (side == 0 && charge > 0) categories.emplace_back("side0_qplus");
    if (side == 0 && charge < 0) categories.emplace_back("side0_qminus");
    if (side == 1 && charge > 0) categories.emplace_back("side1_qplus");
    if (side == 1 && charge < 0) categories.emplace_back("side1_qminus");

    return categories;
  }
}

void MakeResidualHistograms(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputDir = "output",
  const char* outputName = "residual_histograms.root",
  const char* treeName = "residuals",
  const bool recalc_rDCA_zero = true
)
{
  TH1::AddDirectory(kTRUE);

  const TString chainPattern =
    TString::Format("%s/%s", inputDir, filePattern);

  TChain chain(treeName);
  const int nFiles = chain.Add(chainPattern);

  if (nFiles <= 0)
  {
    std::cerr << "ERROR: no files matched: " << chainPattern << std::endl;
    return;
  }

  std::cout << "Added " << nFiles << " files to TChain '" << treeName << "'.\n"
            << "Total entries: " << chain.GetEntries() << std::endl;

  // Scalar branches.
  Int_t side = -1;
  UInt_t ntpc_clusters = 0;

  Double_t pt = 0.;
  Double_t px = 0.;
  Double_t py = 0.;
  Double_t pz = 0.;
  Double_t eta = 0.;
  Double_t charge = 0.;
  Double_t dedx = 0.;
  Double_t vertex_z = 0.;
  Double_t pca_z = 0.;
  Double_t rDCA_zero = 0.;
  Double_t zDCA = 0.;
  Double_t quality = 0.;
  Double_t pca_x = 0.;
  Double_t pca_y = 0.;

  Double_t new_vertex_x = 0.158;
  Double_t new_vertex_y = 0.285;

  // Vector branches.
  std::vector<unsigned int>* layer = nullptr;
  std::vector<double>* cluster_r = nullptr;
  std::vector<double>* cluster_phi = nullptr;
  std::vector<double>* cluster_z = nullptr;
  std::vector<double>* residual_rphi = nullptr;
  std::vector<double>* residual_z = nullptr;

  chain.SetBranchAddress("side", &side);
  chain.SetBranchAddress("ntpc_clusters", &ntpc_clusters);
  chain.SetBranchAddress("pt", &pt);
  chain.SetBranchAddress("px", &px);
  chain.SetBranchAddress("py", &py);
  chain.SetBranchAddress("pz", &pz);
  chain.SetBranchAddress("eta", &eta);
  chain.SetBranchAddress("charge", &charge);
  chain.SetBranchAddress("dedx", &dedx);
  chain.SetBranchAddress("vertex_z", &vertex_z);
  chain.SetBranchAddress("pca_z", &pca_z);
  chain.SetBranchAddress("rDCA_zero", &rDCA_zero);
  chain.SetBranchAddress("zDCA", &zDCA);
  chain.SetBranchAddress("quality", &quality);
  chain.SetBranchAddress("layer", &layer);
  chain.SetBranchAddress("cluster_r", &cluster_r);
  chain.SetBranchAddress("residual_rphi", &residual_rphi);
  chain.SetBranchAddress("residual_z", &residual_z);

  if (recalc_rDCA_zero)
  {
    chain.SetBranchAddress("pca_x", &pca_x);
    chain.SetBranchAddress("pca_y", &pca_y);
  }

  const bool hasClusterPhi = (chain.GetBranch("cluster_phi") != nullptr);
  if (hasClusterPhi)
  {
    chain.SetBranchAddress("cluster_phi", &cluster_phi);
  }
  else
  {
    std::cerr << "WARNING: branch 'cluster_phi' is not available. "
              << "Layer-vs-cluster-phi residual profiles will remain empty."
              << std::endl;
  }

  const bool hasClusterZ = (chain.GetBranch("cluster_z") != nullptr);
  if (hasClusterZ)
  {
    chain.SetBranchAddress("cluster_z", &cluster_z);
  }
  else
  {
    std::cerr << "WARNING: branch 'cluster_z' is not available. "
              << "Cluster-z-vs-layer residual profiles will remain empty."
              << std::endl;
  }

  gSystem->mkdir(outputDir, kTRUE);
  const TString outputPath =
    TString::Format("%s/%s", outputDir, outputName);

  std::unique_ptr<TFile> output(TFile::Open(outputPath, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: could not create output file: "
              << outputPath << std::endl;
    return;
  }

  const std::vector<std::string> categoryNames = {
    "all",
    "side0", "side1",
    "qplus", "qminus",
    "side0_qplus", "side0_qminus",
    "side1_qplus", "side1_qminus"
  };

  std::map<std::string, HistSet> histograms;
  for (const auto& category : categoryNames)
  {
    TDirectory* dir = output->mkdir(category.c_str());
    histograms.emplace(category, bookHistograms(dir, category));
  }

  const Long64_t nEntries = chain.GetEntries();
  for (Long64_t entry = 0; entry < nEntries; ++entry)
  {
    chain.GetEntry(entry);
    
    if (entry % 100000 == 0)
    {
      std::cout << "Processing entry " << entry
                << " / " << nEntries << std::endl;
    }

    if (!std::isfinite(pt) || !std::isfinite(px) || !std::isfinite(py) ||
        !std::isfinite(pz) || !std::isfinite(charge))
    {
      continue;
    }

    // Keep the original baseline requirement for all pre-existing plots.
    // The two new loose residual maps below are allowed to bypass it.

    if (std::abs(pca_z) >= 10. || quality > 20.0)
    {
      continue;
    }

    if (recalc_rDCA_zero)
    {
      rDCA_zero = std::hypot(pca_x + py/(0.003*charge*1.4) - new_vertex_x, pca_y - px/(0.003*charge*1.4) - new_vertex_y) - std::abs(pt/(0.003*charge*1.4));
    }

    const bool passBaseline = ntpc_clusters >= 20;

    const auto activeCategories = categoriesForTrack(side, charge);
    const double phi_px_py = std::atan2(px, py);
    const double qOverPt = (pt != 0.) ? charge / pt : 0.;
    const double p = std::sqrt(pt * pt + pz * pz);
    const double signedP = charge * p;

    // Original selection:
    // ntpc_clusters>30 && pt>0.5 && abs(zDCA)<2
    const bool passSidePca =
      passBaseline && ntpc_clusters > 20 && pt > 0.2 && std::abs(zDCA) < 2e6;

    // Original selection:
    // ntpc_clusters>30 && pca_z in (-15,15) && pt>0.3
    const bool passDcaPhiEta =
      passBaseline &&
      ntpc_clusters > 30 &&
      std::abs(pca_z) < 15. &&
      pt > 0.3;

    // Original q/pT plot had no explicit pT threshold.
    const bool passRDcaQOverPt =
      passBaseline &&
      ntpc_clusters > 30 &&
      std::abs(pca_z) < 15. &&
      pt > 0.;

    // New phase-space QA for studying charge symmetry and distortion stability.
    // These plots use the same fixed PCA-z region on both sides.
    const bool passDistortionPhaseSpaceQA =
      ntpc_clusters > 20 &&
      std::abs(pca_z) < 10. &&
      pt > 0.2 &&
      std::isfinite(eta) &&
      std::isfinite(rDCA_zero) &&
      std::abs(rDCA_zero) < 99.;

    // Original dE/dx selection.
    const bool passDedx =
      passBaseline &&
      ntpc_clusters > 30 &&
      std::abs(rDCA_zero) < 2. &&
      std::abs(zDCA) < 2;

    // All cluster-residual plots require a good transverse DCA.
    // The rDCA_zero distributions themselves remain uncut in rDCA_zero.
    const bool passResidualTrack =
      passBaseline &&
      std::isfinite(rDCA_zero) && std::abs(rDCA_zero) < 10.;

    // New loose residual-map selection only. It intentionally bypasses
    // the baseline ntpc_clusters >= 30 selection used by all existing plots.
    const bool passClusterMapLoose =
      ntpc_clusters > 20 &&
      pt > 0.2 &&
      std::abs(pca_z) < 10. &&
      std::isfinite(rDCA_zero) &&
      std::abs(rDCA_zero) < 2.;

    // Existing cluster residual histograms retain their pT > 2 GeV/c cut.
    const bool passClusterResidual = passResidualTrack && pt > 0.2;

    // New residual maps are made with two thresholds.
    const bool passClusterMapPt1 = passResidualTrack && pt > 0.2;
    const bool passClusterMapPt2 = passResidualTrack && pt > 2.;

    for (const auto& category : activeCategories)
    {
      HistSet& h = histograms.at(category);

      if (passSidePca)
      {
        h.h_side_pca_z->Fill(pca_z, side);
        h.h_side_pca_z_pt->Fill(pca_z, side, pt);
      }

      if (passDcaPhiEta)
      {
        h.h_rDCAzero_phi->Fill(phi_px_py, rDCA_zero);
        h.h_rDCAzero_phi_pt->Fill(phi_px_py, rDCA_zero, pt);

        h.h_zDCA_phi->Fill(phi_px_py, zDCA);
        h.h_zDCA_phi_pt->Fill(phi_px_py, zDCA, pt);

        h.h_zDCA_eta->Fill(eta, zDCA);
        h.h_zDCA_eta_pt->Fill(eta, zDCA, pt);
      }

      if (passRDcaQOverPt)
      {
        h.h_rDCA_qOverPt->Fill(qOverPt, rDCA_zero);
        h.h_rDCA_qOverPt_pt->Fill(qOverPt, rDCA_zero, pt);
      }

      if (passDistortionPhaseSpaceQA)
      {
        h.h_pt_phi_rdca3->Fill(phi_px_py, pt);
        h.h_pt_eta_rdca3->Fill(eta, pt);

        h.h_rDCAzero_phi_rdca3->Fill(phi_px_py, rDCA_zero);
        h.h_rDCAzero_eta_rdca3->Fill(eta, rDCA_zero);

        h.h_rDCAzero_qOverPt_phi_rdca3->Fill(phi_px_py, qOverPt, rDCA_zero);
        h.h_rDCAzero_qOverPt_eta_rdca3->Fill(eta, qOverPt, rDCA_zero);

        h.p_mean_rDCAzero_phi_eta_rdca3->Fill(phi_px_py, eta, rDCA_zero);
        h.p_mean_qOverPt_phi_eta_rdca3->Fill(phi_px_py, eta, qOverPt);

        
        h.h_pt_phi_eta_rdca3->Fill( phi_px_py, eta, pt);
      }

      if (passDedx)
      {
        h.h_dedx_signedP->Fill(signedP, dedx);
        h.h_dedx_signedP_pt->Fill(signedP, dedx, pt);
      }
    }

    if ((!passClusterResidual && !passClusterMapPt1 && !passClusterMapLoose) ||
        !cluster_r || !residual_rphi || !residual_z || !layer)
    {
      continue;
    }

    std::size_t nClusters = std::min(
      {cluster_r->size(), residual_rphi->size(),
       residual_z->size(), layer->size()});

    if (hasClusterPhi && cluster_phi)
    {
      nClusters = std::min(nClusters, cluster_phi->size());
    }
    if (hasClusterZ && cluster_z)
    {
      nClusters = std::min(nClusters, cluster_z->size());
    }

    for (std::size_t i = 0; i < nClusters; ++i)
    {
      const double radius = cluster_r->at(i);
      const double drphi = residual_rphi->at(i);
      const double dz = residual_z->at(i);
      const double tpcLayer = static_cast<double>(layer->at(i));
      const double clusterPhi =
        (hasClusterPhi && cluster_phi) ? cluster_phi->at(i) : 0.;
      const double clusterZ =
        (hasClusterZ && cluster_z) ? cluster_z->at(i) : 0.;

      if (!std::isfinite(radius) || !std::isfinite(drphi) ||
          !std::isfinite(dz) ||
          (hasClusterPhi && !std::isfinite(clusterPhi)) ||
          (hasClusterZ && !std::isfinite(clusterZ)))
      {
        continue;
      }

      // Reject pathological residual outliers at the filling stage.
      // Apply the rejection independently so a bad rphi residual does not
      // discard an otherwise valid z residual, and vice versa.
      const bool goodRphiResidual = std::abs(drphi) < 2.;
      const bool goodZResidual = std::abs(dz) < 9999.;

      for (const auto& category : activeCategories)
      {
        HistSet& h = histograms.at(category);

        if (passClusterResidual)
        {
          if (goodRphiResidual)
          {
            h.h_residual_rphi_cluster_r->Fill(radius, drphi);
            h.h_residual_rphi_cluster_r_pt->Fill(radius, drphi, pt);
            h.h_residual_rphi_layer->Fill(tpcLayer, drphi);
            h.h_residual_rphi_layer_pt->Fill(tpcLayer, drphi, pt);
          }

          if (goodZResidual)
          {
            h.h_residual_z_cluster_r->Fill(radius, dz);
            h.h_residual_z_cluster_r_pt->Fill(radius, dz, pt);
            h.h_residual_z_layer->Fill(tpcLayer, dz);
            h.h_residual_z_layer_pt->Fill(tpcLayer, dz, pt);
          }
        }

        if (hasClusterPhi && cluster_phi && passClusterMapLoose)
        {
          if (goodRphiResidual)
          {
            h.h_mean_residual_rphi_layer_phi_pt0->Fill(tpcLayer, clusterPhi, drphi);
          }
          if (goodZResidual)
          {
            h.h_mean_residual_z_layer_phi_pt0->Fill(tpcLayer, clusterPhi, dz);
          }
        }

        if (hasClusterPhi && cluster_phi && passClusterMapPt1)
        {
          if (goodRphiResidual)
          {
            h.h_mean_residual_rphi_layer_phi_pt1->Fill(tpcLayer, clusterPhi, drphi);
          }
          if (goodZResidual)
          {
            h.h_mean_residual_z_layer_phi_pt1->Fill(tpcLayer, clusterPhi, dz);
          }
        }

        if (hasClusterPhi && cluster_phi && passClusterMapPt2)
        {
          if (goodRphiResidual)
          {
            h.h_mean_residual_rphi_layer_phi_pt2->Fill(tpcLayer, clusterPhi, drphi);
          }
          if (goodZResidual)
          {
            h.h_mean_residual_z_layer_phi_pt2->Fill(tpcLayer, clusterPhi, dz);
          }
        }

        if (hasClusterZ && cluster_z && goodZResidual && passClusterMapPt1)
        {
          h.h_mean_residual_z_cluster_z_layer_pt1->Fill(clusterZ, tpcLayer, dz);
        }

        if (hasClusterZ && cluster_z && goodZResidual && passClusterMapPt2)
        {
          h.h_mean_residual_z_cluster_z_layer_pt2->Fill(clusterZ, tpcLayer, dz);
        }
      }
    }
  }

  output->cd();
  output->Write();
  output->Close();

  std::cout << "Wrote histograms to: " << outputPath << std::endl;
}