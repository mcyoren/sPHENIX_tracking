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

    // Mean residual maps versus detector layer and cluster phi.
    // Two pT selections are stored separately: pT > 1 and pT > 2 GeV/c.
    // Option "s" makes profile-bin errors equal to the residual spread
    // (standard deviation), rather than the error on the mean.
    TProfile2D* h_mean_residual_rphi_layer_phi_pt1 = nullptr;
    TProfile2D* h_mean_residual_rphi_layer_phi_pt2 = nullptr;
    TProfile2D* h_mean_residual_z_layer_phi_pt1 = nullptr;
    TProfile2D* h_mean_residual_z_layer_phi_pt2 = nullptr;
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

    h.h_dedx_signedP =
      new TH2D("h_dedx_vs_signedP",
               "dE/dx vs charge#times(1-2side)#timesp" + suffix +
               ";q(1-2side)p [GeV/c];dE/dx",
               400, -2., 2., 1000, 0., 6000.);
    h.h_dedx_signedP_pt =
      new TH3D("h_dedx_vs_signedP_vs_pt",
               "dE/dx vs signed p vs p_{T}" + suffix +
               ";q(1-2side)p [GeV/c];dE/dx;p_{T} [GeV/c]",
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
  const char* treeName = "residuals")
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

  // Vector branches.
  std::vector<unsigned int>* layer = nullptr;
  std::vector<double>* cluster_r = nullptr;
  std::vector<double>* cluster_phi = nullptr;
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

  chain.SetBranchAddress("layer", &layer);
  chain.SetBranchAddress("cluster_r", &cluster_r);
  chain.SetBranchAddress("residual_rphi", &residual_rphi);
  chain.SetBranchAddress("residual_z", &residual_z);

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

    // Global track-quality selection applied to every histogram below.
    if (std::abs(pca_z) > 10. ||
        ntpc_clusters < 30 ||
        !std::isfinite(rDCA_zero) ||
        std::abs(rDCA_zero) >= 2.)
    {
      continue;
    }

    const auto activeCategories = categoriesForTrack(side, charge);
    const double phi_px_py = std::atan2(px, py);
    const double qOverPt = (pt != 0.) ? charge / pt : 0.;
    const double p = std::sqrt(pt * pt + pz * pz);
    const double signedP = charge * (1. - 0. * side) * p;

    // Original selection:
    // ntpc_clusters>30 && pt>0.5 && abs(zDCA)<2
    const bool passSidePca =
      ntpc_clusters > 30 && pt > 0.5 && std::abs(zDCA) < 2e6;

    // Original selection:
    // ntpc_clusters>30 && pca_z in (-15,15) && pt>0.3
    const bool passDcaPhiEta =
      ntpc_clusters > 30 &&
      std::abs(pca_z) < 15. &&
      pt > 0.3;

    // Original q/pT plot had no explicit pT threshold.
    const bool passRDcaQOverPt =
      ntpc_clusters > 30 &&
      std::abs(pca_z) < 15. &&
      pt > 0.;

    // Original dE/dx selection.
    const bool passDedx =
      ntpc_clusters > 25 &&
      std::abs(rDCA_zero) < 2. &&
      std::abs(zDCA) < 2;

    // Existing cluster residual histograms retain their pT > 2 GeV/c cut.
    const bool passClusterResidual = pt > 2.;

    // New layer-vs-cluster-phi maps are made with two thresholds.
    const bool passClusterMapPt1 = pt > 1.;
    const bool passClusterMapPt2 = pt > 2.;

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

      if (passDedx)
      {
        h.h_dedx_signedP->Fill(signedP, dedx);
        h.h_dedx_signedP_pt->Fill(signedP, dedx, pt);
      }
    }

    if ((!passClusterResidual && !passClusterMapPt1) ||
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

    for (std::size_t i = 0; i < nClusters; ++i)
    {
      const double radius = cluster_r->at(i);
      const double drphi = residual_rphi->at(i);
      const double dz = residual_z->at(i);
      const double tpcLayer = static_cast<double>(layer->at(i));
      const double clusterPhi =
        (hasClusterPhi && cluster_phi) ? cluster_phi->at(i) : 0.;

      if (!std::isfinite(radius) || !std::isfinite(drphi) ||
          !std::isfinite(dz) ||
          (hasClusterPhi && !std::isfinite(clusterPhi)))
      {
        continue;
      }

      for (const auto& category : activeCategories)
      {
        HistSet& h = histograms.at(category);

        if (passClusterResidual)
        {
          h.h_residual_rphi_cluster_r->Fill(radius, drphi);
          h.h_residual_rphi_cluster_r_pt->Fill(radius, drphi, pt);

          h.h_residual_z_cluster_r->Fill(radius, dz);
          h.h_residual_z_cluster_r_pt->Fill(radius, dz, pt);

          h.h_residual_rphi_layer->Fill(tpcLayer, drphi);
          h.h_residual_rphi_layer_pt->Fill(tpcLayer, drphi, pt);

          h.h_residual_z_layer->Fill(tpcLayer, dz);
          h.h_residual_z_layer_pt->Fill(tpcLayer, dz, pt);
        }

        if (hasClusterPhi && cluster_phi && passClusterMapPt1)
        {
          h.h_mean_residual_rphi_layer_phi_pt1->Fill(tpcLayer, clusterPhi, drphi);
          h.h_mean_residual_z_layer_phi_pt1->Fill(tpcLayer, clusterPhi, dz);
        }

        if (hasClusterPhi && cluster_phi && passClusterMapPt2)
        {
          h.h_mean_residual_rphi_layer_phi_pt2->Fill(tpcLayer, clusterPhi, drphi);
          h.h_mean_residual_z_layer_phi_pt2->Fill(tpcLayer, clusterPhi, dz);
        }
      }
    }
  }

  output->cd();
  output->Write();
  output->Close();

  std::cout << "Wrote histograms to: " << outputPath << std::endl;
}
