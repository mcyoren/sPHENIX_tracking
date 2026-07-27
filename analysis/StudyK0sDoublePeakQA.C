// StudyK0sDoublePeakQA.C
//
// One-pass QA macro for diagnosing a double structure in K0S -> pi+ pi-.
// It applies exactCut1 from MakeMineK0sPairHistograms.C and fills all QA
// histograms during the same TChain loop.
//
// Example:
// root -l -b -q 'StudyK0sDoublePeakQA.C("/path/to/files","pair*.root","output","k0s_double_peak_qa.root","pairTree")'
//
// Optional peak-window arguments can be changed at the end of the call:
//   low peak  = [lowMassMin,  lowMassMax]
//   high peak = [highMassMin, highMassMax]

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
  constexpr double kPionMass = 0.13957039;
  constexpr double kTwoPi = 2.0 * TMath::Pi();

  double wrapPhi(double phi)
  {
    while (phi >= TMath::Pi()) phi -= kTwoPi;
    while (phi < -TMath::Pi()) phi += kTwoPi;
    return phi;
  }

  double safeEta(double px, double py, double pz)
  {
    const double pt = std::hypot(px, py);
    return (pt > 0.) ? std::asinh(pz / pt) : 0.;
  }

  double safeMass(double px1, double py1, double pz1,
                  double px2, double py2, double pz2)
  {
    const double p1sq = px1*px1 + py1*py1 + pz1*pz1;
    const double p2sq = px2*px2 + py2*py2 + pz2*pz2;
    const double e1 = std::sqrt(p1sq + kPionMass*kPionMass);
    const double e2 = std::sqrt(p2sq + kPionMass*kPionMass);
    const double px = px1 + px2;
    const double py = py1 + py2;
    const double pz = pz1 + pz2;
    const double m2 = (e1 + e2)*(e1 + e2) - px*px - py*py - pz*pz;
    return (m2 > 0.) ? std::sqrt(m2) : 0.;
  }

  double openingAngle(double px1, double py1, double pz1,
                      double px2, double py2, double pz2)
  {
    const double p1 = std::sqrt(px1*px1 + py1*py1 + pz1*pz1);
    const double p2 = std::sqrt(px2*px2 + py2*py2 + pz2*pz2);
    if (p1 <= 0. || p2 <= 0.) return -1.;
    double c = (px1*px2 + py1*py2 + pz1*pz2)/(p1*p2);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c);
  }

  bool branchExists(TChain& chain, const char* name)
  {
    return chain.GetBranch(name) != nullptr;
  }

  TH2F* bookMassVs(TDirectory* dir,
                   const char* name,
                   const char* xTitle,
                   int nx, double xmin, double xmax,
                   int nmass = 300, double mmin = 0.42, double mmax = 0.57)
  {
    dir->cd();
    return new TH2F(name,
                    TString::Format("K^{0}_{S} mass vs %s; %s; m_{#pi#pi} [GeV/c^{2}]",
                                    xTitle, xTitle),
                    nx, xmin, xmax, nmass, mmin, mmax);
  }

  struct MassVsQA
  {
    TH2F* v0Pt = nullptr;
    TH2F* v0Phi = nullptr;
    TH2F* v0Eta = nullptr;
    TH2F* decayR = nullptr;
    TH2F* decayPhi = nullptr;
    TH2F* pcaZ = nullptr;
    TH2F* deltaPcaZ = nullptr;
    TH2F* pairDCA = nullptr;
    TH2F* dira = nullptr;
    TH2F* alpha = nullptr;
    TH2F* absAlpha = nullptr;

    TH2F* ptPos = nullptr;
    TH2F* ptNeg = nullptr;
    TH2F* ptMin = nullptr;
    TH2F* ptMax = nullptr;
    TH2F* ptAsymmetry = nullptr;
    TH2F* etaPos = nullptr;
    TH2F* etaNeg = nullptr;
    TH2F* etaMin = nullptr;
    TH2F* etaMax = nullptr;
    TH2F* phiPos = nullptr;
    TH2F* phiNeg = nullptr;
    TH2F* deltaPhi = nullptr;
    TH2F* opening = nullptr;
    TH2F* qOverPtPos = nullptr;
    TH2F* qOverPtNeg = nullptr;
    TH2F* pzPos = nullptr;
    TH2F* pzNeg = nullptr;
    TH2F* npointsMin = nullptr;
    TH2F* npointsMax = nullptr;
    TH2F* qualityMin = nullptr;
    TH2F* qualityMax = nullptr;
    TH2F* massRecalc = nullptr;
    TH2F* massDifference = nullptr;
    TH3F* deltaPhiVsMassVsPt = nullptr;
    TH2F* massvsPtminusoverPtplus = nullptr;
  };

  MassVsQA bookMassVsQA(TDirectory* dir)
  {
    MassVsQA h;
    h.v0Pt          = bookMassVs(dir, "h_mass_vs_v0_pt", "p_{T}^{V0} [GeV/c]", 100, 0., 5.);
    h.v0Phi         = bookMassVs(dir, "h_mass_vs_v0_phi", "#phi_{V0} [rad]", 144, -TMath::Pi(), TMath::Pi());
    h.v0Eta         = bookMassVs(dir, "h_mass_vs_v0_eta", "#eta_{V0}", 120, -1.5, 1.5);
    h.decayR        = bookMassVs(dir, "h_mass_vs_decay_radius", "decay radius [cm]", 120, 0., 80.);
    h.decayPhi      = bookMassVs(dir, "h_mass_vs_decay_phi", "#phi_{decay} [rad]", 144, -TMath::Pi(), TMath::Pi());
    h.pcaZ          = bookMassVs(dir, "h_mass_vs_pca_z", "PCA_{z} [cm]", 120, -15., 15.);
    h.deltaPcaZ     = bookMassVs(dir, "h_mass_vs_abs_delta_pca_z", "|PCA_{z,1}-PCA_{z,2}| [cm]", 100, 0., 0.5);
    h.pairDCA       = bookMassVs(dir, "h_mass_vs_abs_pairDCA", "|pair DCA| [cm]", 100, 0., 1.);
    h.dira          = bookMassVs(dir, "h_mass_vs_DIRA", "DIRA", 100, 0.8, 1.0);
    h.alpha         = bookMassVs(dir, "h_mass_vs_alpha", "#alpha", 120, -0.8, 0.8);
    h.absAlpha      = bookMassVs(dir, "h_mass_vs_abs_alpha", "|#alpha|", 80, 0., 0.8);

    h.ptPos         = bookMassVs(dir, "h_mass_vs_pt_positive", "p_{T}^{+} [GeV/c]", 100, 0., 5.);
    h.ptNeg         = bookMassVs(dir, "h_mass_vs_pt_negative", "p_{T}^{-} [GeV/c]", 100, 0., 5.);
    h.ptMin         = bookMassVs(dir, "h_mass_vs_pt_min", "min(p_{T}^{+},p_{T}^{-}) [GeV/c]", 100, 0., 3.);
    h.ptMax         = bookMassVs(dir, "h_mass_vs_pt_max", "max(p_{T}^{+},p_{T}^{-}) [GeV/c]", 100, 0., 5.);
    h.ptAsymmetry   = bookMassVs(dir, "h_mass_vs_pt_asymmetry", "(p_{T}^{+}-p_{T}^{-})/(p_{T}^{+}+p_{T}^{-})", 120, -1., 1.);
    h.etaPos        = bookMassVs(dir, "h_mass_vs_eta_positive", "#eta^{+}", 120, -1.5, 1.5);
    h.etaNeg        = bookMassVs(dir, "h_mass_vs_eta_negative", "#eta^{-}", 120, -1.5, 1.5);
    h.etaMin        = bookMassVs(dir, "h_mass_vs_eta_min", "min(#eta^{+},#eta^{-})", 120, -1.5, 1.5);
    h.etaMax        = bookMassVs(dir, "h_mass_vs_eta_max", "max(#eta^{+},#eta^{-})", 120, -1.5, 1.5);
    h.phiPos        = bookMassVs(dir, "h_mass_vs_phi_positive", "#phi^{+} [rad]", 144, -TMath::Pi(), TMath::Pi());
    h.phiNeg        = bookMassVs(dir, "h_mass_vs_phi_negative", "#phi^{-} [rad]", 144, -TMath::Pi(), TMath::Pi());
    h.deltaPhi      = bookMassVs(dir, "h_mass_vs_delta_phi", "#Delta#phi_{+-} [rad]", 144, -TMath::Pi(), TMath::Pi());
    h.opening       = bookMassVs(dir, "h_mass_vs_opening_angle", "opening angle [rad]", 120, 0., TMath::Pi());
    h.qOverPtPos    = bookMassVs(dir, "h_mass_vs_qOverPt_positive", "q/p_{T}^{+} [(GeV/c)^{-1}]", 120, 0., 4.);
    h.qOverPtNeg    = bookMassVs(dir, "h_mass_vs_qOverPt_negative", "q/p_{T}^{-} [(GeV/c)^{-1}]", 120, -4., 0.);
    h.pzPos         = bookMassVs(dir, "h_mass_vs_pz_positive", "p_{z}^{+} [GeV/c]", 120, -5., 5.);
    h.pzNeg         = bookMassVs(dir, "h_mass_vs_pz_negative", "p_{z}^{-} [GeV/c]", 120, -5., 5.);
    h.npointsMin    = bookMassVs(dir, "h_mass_vs_npoints_min", "min(N_{points}^{+},N_{points}^{-})", 40, 30., 70.);
    h.npointsMax    = bookMassVs(dir, "h_mass_vs_npoints_max", "max(N_{points}^{+},N_{points}^{-})", 40, 30., 70.);
    h.qualityMin    = bookMassVs(dir, "h_mass_vs_quality_min", "min(quality^{+},quality^{-})", 100, 0., 2.);
    h.qualityMax    = bookMassVs(dir, "h_mass_vs_quality_max", "max(quality^{+},quality^{-})", 100, 0., 2.);
    h.massRecalc    = bookMassVs(dir, "h_stored_mass_vs_recalculated_mass", "recalculated m_{#pi#pi} [GeV/c^{2}]", 300, 0.42, 0.57);
    h.deltaPhiVsMassVsPt = new TH3F("h_delta_phi_vs_mass_vs_pt", "K^{0}_{S} #Delta#phi_{+-} vs m_{#pi#pi} vs p_{T}^{V0}; m_{#pi#pi} [GeV/c^{2}]; p_{T}^{V0} [GeV/c]; #Delta#phi_{+-} [rad]",
      300, 0.42, 0.57, 50, 0., 5., 144, -TMath::Pi(), TMath::Pi());
    h.massDifference= bookMassVs(dir, "h_mass_vs_stored_minus_recalculated", "stored mass - recalculated mass [GeV/c^{2}]", 200, -0.02, 0.02);
    h.massvsPtminusoverPtplus = bookMassVs(dir, "h_mass_vs_pt_minus_over_pt_plus", "p_{T}^{-}/p_{T}^{+}", 100, 0., 5.0);
    return h;
  }

  void fillMassVsQA(MassVsQA& h,
                    double mass,
                    double v0Pt, double v0Phi, double v0Eta,
                    double decayR, double decayPhi, double pcaZ,
                    double deltaPcaZ, double pairDCA, double dira, double alpha,
                    double ptPos, double ptNeg,
                    double etaPos, double etaNeg,
                    double phiPos, double phiNeg,
                    double pzPos, double pzNeg,
                    int npointsPos, int npointsNeg,
                    double qualityPos, double qualityNeg,
                    double opening, double massRecalc)
  {
    const double ptSum = ptPos + ptNeg;
    const double ptAsym = (ptSum > 0.) ? (ptPos - ptNeg)/ptSum : 0.;

    h.v0Pt->Fill(v0Pt, mass);
    h.v0Phi->Fill(v0Phi, mass);
    h.v0Eta->Fill(v0Eta, mass);
    h.decayR->Fill(decayR, mass);
    h.decayPhi->Fill(decayPhi, mass);
    h.pcaZ->Fill(pcaZ, mass);
    h.deltaPcaZ->Fill(deltaPcaZ, mass);
    h.pairDCA->Fill(pairDCA, mass);
    h.dira->Fill(dira, mass);
    h.alpha->Fill(alpha, mass);
    h.absAlpha->Fill(std::abs(alpha), mass);

    h.ptPos->Fill(ptPos, mass);
    h.ptNeg->Fill(ptNeg, mass);
    h.ptMin->Fill(std::min(ptPos, ptNeg), mass);
    h.ptMax->Fill(std::max(ptPos, ptNeg), mass);
    h.ptAsymmetry->Fill(ptAsym, mass);
    h.etaPos->Fill(etaPos, mass);
    h.etaNeg->Fill(etaNeg, mass);
    h.etaMin->Fill(std::min(etaPos, etaNeg), mass);
    h.etaMax->Fill(std::max(etaPos, etaNeg), mass);
    h.phiPos->Fill(phiPos, mass);
    h.phiNeg->Fill(phiNeg, mass);
    ///h.deltaPhi->Fill(wrapPhi(phiPos - phiNeg), mass);
    h.opening->Fill(opening, mass);
    if (ptPos > 0.) h.qOverPtPos->Fill(1.0/ptPos, mass);
    if (ptNeg > 0.) h.qOverPtNeg->Fill(-1.0/ptNeg, mass);
    h.pzPos->Fill(pzPos, mass);
    h.pzNeg->Fill(pzNeg, mass);
    h.npointsMin->Fill(std::min(npointsPos, npointsNeg), mass);
    h.npointsMax->Fill(std::max(npointsPos, npointsNeg), mass);
    h.qualityMin->Fill(std::min(qualityPos, qualityNeg), mass);
    h.qualityMax->Fill(std::max(qualityPos, qualityNeg), mass);
    h.massRecalc->Fill(massRecalc, mass);
    h.massDifference->Fill(mass - massRecalc, mass);
    h.massvsPtminusoverPtplus->Fill(ptNeg / ptPos, mass);
    //h.deltaPhiVsMassVsPt->Fill(mass, v0Pt, wrapPhi(phiPos - phiNeg));
  }

  struct WindowMaps
  {
    TH2F* phiPosVsPhiNeg = nullptr;
    TH2F* ptPosVsPtNeg = nullptr;
    TH2F* etaPosVsEtaNeg = nullptr;
    TH2F* decayPhiVsR = nullptr;
    TH2F* v0PhiVsPt = nullptr;
    TH2F* phiPosVsPtPos = nullptr;
    TH2F* phiNegVsPtNeg = nullptr;
  };

  WindowMaps bookWindowMaps(TDirectory* dir, const char* label)
  {
    dir->cd();
    WindowMaps h;
    h.phiPosVsPhiNeg = new TH2F("h_phi_positive_vs_phi_negative",
      TString::Format("%s;#phi^{-} [rad];#phi^{+} [rad]", label),
      144, -TMath::Pi(), TMath::Pi(), 144, -TMath::Pi(), TMath::Pi());
    h.ptPosVsPtNeg = new TH2F("h_pt_positive_vs_pt_negative",
      TString::Format("%s;p_{T}^{-} [GeV/c];p_{T}^{+} [GeV/c]", label),
      100, 0., 4., 100, 0., 4.);
    h.etaPosVsEtaNeg = new TH2F("h_eta_positive_vs_eta_negative",
      TString::Format("%s;#eta^{-};#eta^{+}", label),
      120, -1.5, 1.5, 120, -1.5, 1.5);
    h.decayPhiVsR = new TH2F("h_decay_phi_vs_radius",
      TString::Format("%s;decay radius [cm];#phi_{decay} [rad]", label),
      120, 0., 80., 144, -TMath::Pi(), TMath::Pi());
    h.v0PhiVsPt = new TH2F("h_v0_phi_vs_v0_pt",
      TString::Format("%s;p_{T}^{V0} [GeV/c];#phi_{V0} [rad]", label),
      100, 0., 5., 144, -TMath::Pi(), TMath::Pi());
    h.phiPosVsPtPos = new TH2F("h_phi_positive_vs_pt_positive",
      TString::Format("%s;p_{T}^{+} [GeV/c];#phi^{+} [rad]", label),
      100, 0., 4., 144, -TMath::Pi(), TMath::Pi());
    h.phiNegVsPtNeg = new TH2F("h_phi_negative_vs_pt_negative",
      TString::Format("%s;p_{T}^{-} [GeV/c];#phi^{-} [rad]", label),
      100, 0., 4., 144, -TMath::Pi(), TMath::Pi());
    return h;
  }

  void fillWindowMaps(WindowMaps& h,
                      double phiPos, double phiNeg,
                      double ptPos, double ptNeg,
                      double etaPos, double etaNeg,
                      double decayR, double decayPhi,
                      double v0Pt, double v0Phi)
  {
    h.phiPosVsPhiNeg->Fill(phiNeg, phiPos);
    h.ptPosVsPtNeg->Fill(ptNeg, ptPos);
    h.etaPosVsEtaNeg->Fill(etaNeg, etaPos);
    h.decayPhiVsR->Fill(decayR, decayPhi);
    h.v0PhiVsPt->Fill(v0Pt, v0Phi);
    h.phiPosVsPtPos->Fill(ptPos, phiPos);
    h.phiNegVsPtNeg->Fill(ptNeg, phiNeg);
  }
}

void StudyK0sDoublePeakQA(const char* inputDir = ".",
                          const char* filePattern = "*.root",
                          const char* outputDir = "output",
                          const char* outputName = "k0s_double_peak_qa.root",
                          const char* treeName = "pairTree",
                          double lowMassMin = 0.470,
                          double lowMassMax = 0.490,
                          double highMassMin = 0.490,
                          double highMassMax = 0.510,
                          Long64_t maxEntries = -1)
{
  TH1::AddDirectory(kTRUE);

  const TString chainPattern = TString::Format("%s/%s", inputDir, filePattern);
  TChain chain(treeName);
  const int nFiles = chain.Add(chainPattern);
  if (nFiles <= 0)
  {
    std::cerr << "ERROR: no files matched " << chainPattern << std::endl;
    return;
  }

  const std::vector<std::string> requiredBranches = {
    "mass_Kshort", "v0_pt", "pca_x", "pca_y", "pca_z",
    "pca1_z", "pca2_z", "px1", "py1", "pz1", "px2", "py2", "pz2",
    "v0_px", "v0_py", "v0_pz", "alpha", "pairDCA",
    "npoints1", "npoints2", "charge1", "charge2", "quality1", "quality2"
  };

  bool missing = false;
  for (const auto& name : requiredBranches)
  {
    if (!branchExists(chain, name.c_str()))
    {
      std::cerr << "ERROR: missing branch " << name << std::endl;
      missing = true;
    }
  }
  if (missing) return;

  Float_t mass_Kshort = 0.f, v0_pt = 0.f;
  Float_t pca_x = 0.f, pca_y = 0.f, pca_z = 0.f, pca1_z = 0.f, pca2_z = 0.f;
  Float_t px1 = 0.f, py1 = 0.f, pz1 = 0.f, px2 = 0.f, py2 = 0.f, pz2 = 0.f;
  Float_t v0_px = 0.f, v0_py = 0.f, v0_pz = 0.f;
  Float_t alpha = 0.f, pairDCA = 0.f;
  Float_t charge1 = 0.f, charge2 = 0.f, quality1 = 0.f, quality2 = 0.f;
  Short_t npoints1 = 0, npoints2 = 0;

  chain.SetBranchAddress("mass_Kshort", &mass_Kshort);
  chain.SetBranchAddress("v0_pt", &v0_pt);
  chain.SetBranchAddress("pca_x", &pca_x);
  chain.SetBranchAddress("pca_y", &pca_y);
  chain.SetBranchAddress("pca_z", &pca_z);
  chain.SetBranchAddress("pca1_z", &pca1_z);
  chain.SetBranchAddress("pca2_z", &pca2_z);
  chain.SetBranchAddress("px1", &px1);
  chain.SetBranchAddress("py1", &py1);
  chain.SetBranchAddress("pz1", &pz1);
  chain.SetBranchAddress("px2", &px2);
  chain.SetBranchAddress("py2", &py2);
  chain.SetBranchAddress("pz2", &pz2);
  chain.SetBranchAddress("v0_px", &v0_px);
  chain.SetBranchAddress("v0_py", &v0_py);
  chain.SetBranchAddress("v0_pz", &v0_pz);
  chain.SetBranchAddress("alpha", &alpha);
  chain.SetBranchAddress("pairDCA", &pairDCA);
  chain.SetBranchAddress("charge1", &charge1);
  chain.SetBranchAddress("charge2", &charge2);
  chain.SetBranchAddress("quality1", &quality1);
  chain.SetBranchAddress("quality2", &quality2);
  chain.SetBranchAddress("npoints1", &npoints1);
  chain.SetBranchAddress("npoints2", &npoints2);

  gSystem->mkdir(outputDir, kTRUE);
  const TString outputPath = TString::Format("%s/%s", outputDir, outputName);
  std::unique_ptr<TFile> output(TFile::Open(outputPath, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: cannot create " << outputPath << std::endl;
    return;
  }

  TDirectory* massDir = output->mkdir("mass_vs_phase_space_exactCut1_unlike");
  MassVsQA massQA = bookMassVsQA(massDir);

  output->cd();
  TH1F* hMassUnlike = new TH1F("h_mass_exactCut1_unlike",
    "exactCut1 unlike sign;m_{#pi#pi} [GeV/c^{2}];candidates", 600, 0.42, 0.57);
  TH1F* hMassLike = new TH1F("h_mass_exactCut1_like",
    "exactCut1 like sign;m_{#pi#pi} [GeV/c^{2}];candidates", 600, 0.42, 0.57);
  TH2F* hMassVsV0PtLike = new TH2F("h_mass_vs_v0_pt_exactCut1_like",
    "exactCut1 like sign;p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]",
    100, 0., 5., 300, 0.42, 0.57);

  TDirectory* lowDir = output->mkdir("low_mass_window");
  WindowMaps lowMaps = bookWindowMaps(lowDir,
    TString::Format("low mass %.3f-%.3f GeV/c^{2}", lowMassMin, lowMassMax));

  TDirectory* highDir = output->mkdir("high_mass_window");
  WindowMaps highMaps = bookWindowMaps(highDir,
    TString::Format("high mass %.3f-%.3f GeV/c^{2}", highMassMin, highMassMax));

  const Long64_t totalEntries = chain.GetEntries();
  const Long64_t nEntries = (maxEntries >= 0) ? std::min(totalEntries, maxEntries) : totalEntries;
  std::cout << "Added " << nFiles << " files, processing " << nEntries
            << " / " << totalEntries << " entries" << std::endl;
  std::cout << "Low window:  [" << lowMassMin << ", " << lowMassMax << "]\n"
            << "High window: [" << highMassMin << ", " << highMassMax << "]" << std::endl;

  Long64_t nExactUnlike = 0;
  Long64_t nExactLike = 0;

  for (Long64_t entry = 0; entry < nEntries; ++entry)
  {
    chain.GetEntry(entry);
    if (entry % 1000000 == 0)
      std::cout << "Processing " << entry << " / " << nEntries << std::endl;

    const double pt1 = std::hypot(px1, py1);
    const double pt2 = std::hypot(px2, py2);
    const double absDeltaPcaZ = std::abs(pca1_z - pca2_z);
    const double decayRadius = std::hypot(pca_x, pca_y);
    const double absAlpha = std::abs(alpha);
    const double absPairDCA = std::abs(pairDCA);

    const double pMag = std::sqrt(v0_px*v0_px + v0_py*v0_py + v0_pz*v0_pz);
    const double rMag = std::sqrt(pca_x*pca_x + pca_y*pca_y + pca_z*pca_z);
    const double dira = (pMag > 0. && rMag > 0.)
      ? (v0_px*pca_x + v0_py*pca_y + v0_pz*pca_z)/(pMag*rMag)
      : -2.;

    // exactCut1 copied exactly from MakeMineK0sPairHistograms.C.
    // Note: the current code uses quality1/2 < 20.0.
    const bool exactCut1 =
      pca_z > -15.f && pca_z < 15.f &&
      absDeltaPcaZ < 0.5 &&
      pt1 > 0.2 && pt2 > 0.2 &&
      decayRadius > 2.0 &&
      absAlpha < 0.99 &&
      absPairDCA < 1.5 &&
      dira > 0.88 &&
      npoints1 > 20 && npoints2 > 20 &&
      quality1 < 1.5 && quality2 < 1.5;


    const bool exactCut2 =
      pca_z > -10.f && pca_z < 10.f &&
      absDeltaPcaZ < 0.2 &&
      pt1 > 0.2 && pt2 > 0.2 &&
      decayRadius > 2.0 &&
      absAlpha < 0.99 &&
      absPairDCA < 0.5 &&
      dira > 0.95 &&
      npoints1 > 30 && npoints2 > 30 &&
      quality1 < 0.5 && quality2 < 0.5;


    if (!exactCut1) continue;

    const double phiPos0 = std::atan2(charge1==1 ? py1 : py2, charge1==1 ? px1 : px2);
    const double phiNeg0 = std::atan2(charge1==1 ? py2 : py1, charge1==1 ? px2 : px1);
    const double deltaPhi0 = wrapPhi(phiPos0 - phiNeg0);

    if(exactCut2) 
    {
      massQA.deltaPhi->Fill(deltaPhi0, mass_Kshort);
      massQA.deltaPhiVsMassVsPt->Fill(mass_Kshort, v0_pt, deltaPhi0);
    }

    if (deltaPhi0 <0.8 - 0.4 * ( v0_pt  < 2.0 ? v0_pt : 2.0 ) )
    {
      continue;
    }

    const bool unlike = charge1 * charge2 < 0.f;
    const bool like = charge1 * charge2 > 0.f;
    if (!unlike && !like) continue;

    if (like)
    {
      ++nExactLike;
      hMassLike->Fill(mass_Kshort);
      hMassVsV0PtLike->Fill(v0_pt, mass_Kshort);
      continue;
    }

    ++nExactUnlike;
    hMassUnlike->Fill(mass_Kshort);

    const bool track1Positive = charge1 > 0.f;
    const double pxPos = track1Positive ? px1 : px2;
    const double pyPos = track1Positive ? py1 : py2;
    const double pzPos = track1Positive ? pz1 : pz2;
    const double pxNeg = track1Positive ? px2 : px1;
    const double pyNeg = track1Positive ? py2 : py1;
    const double pzNeg = track1Positive ? pz2 : pz1;
    const double ptPos = std::hypot(pxPos, pyPos);
    const double ptNeg = std::hypot(pxNeg, pyNeg);
    const double phiPos = std::atan2(pyPos, pxPos);
    const double phiNeg = std::atan2(pyNeg, pxNeg);
    const double etaPos = safeEta(pxPos, pyPos, pzPos);
    const double etaNeg = safeEta(pxNeg, pyNeg, pzNeg);
    const int npointsPos = track1Positive ? npoints1 : npoints2;
    const int npointsNeg = track1Positive ? npoints2 : npoints1;
    const double qualityPos = track1Positive ? quality1 : quality2;
    const double qualityNeg = track1Positive ? quality2 : quality1;

    const double v0Phi = std::atan2(v0_py, v0_px);
    const double v0Eta = safeEta(v0_px, v0_py, v0_pz);
    const double decayPhi = std::atan2(pca_y, pca_x);
    const double open = openingAngle(px1, py1, pz1, px2, py2, pz2);
    const double recalculatedMass = safeMass(px1, py1, pz1, px2, py2, pz2);

    fillMassVsQA(massQA,
                 mass_Kshort,
                 v0_pt, v0Phi, v0Eta,
                 decayRadius, decayPhi, pca_z,
                 absDeltaPcaZ, absPairDCA, dira, alpha,
                 ptPos, ptNeg,
                 etaPos, etaNeg,
                 phiPos, phiNeg,
                 pzPos, pzNeg,
                 npointsPos, npointsNeg,
                 qualityPos, qualityNeg,
                 open, recalculatedMass);

    if (mass_Kshort >= lowMassMin && mass_Kshort < lowMassMax)
      fillWindowMaps(lowMaps, phiPos, phiNeg, ptPos, ptNeg, etaPos, etaNeg,
                     decayRadius, decayPhi, v0_pt, v0Phi);

    if (mass_Kshort >= highMassMin && mass_Kshort < highMassMax)
      fillWindowMaps(highMaps, phiPos, phiNeg, ptPos, ptNeg, etaPos, etaNeg,
                     decayRadius, decayPhi, v0_pt, v0Phi);
  }

  output->Write();
  output->Close();

  std::cout << "exactCut1 unlike candidates: " << nExactUnlike << std::endl;
  std::cout << "exactCut1 like candidates:   " << nExactLike << std::endl;
  std::cout << "Wrote output: " << outputPath << std::endl;
}
