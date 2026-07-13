// MakeK0sPairHistograms.C
//
// Example:
// root -l -b -q 'MakeK0sPairHistograms.C("/path/to/files","pair*.root","output","k0s_qa.root","pairTree")'
//
// The TChain is read exactly once.

#include <TBranch.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH3F.h>
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
  enum class ChargeCategory { Unlike, Like, PlusPlus, MinusMinus };

  struct Selection
  {
    std::string name;
    double maxAbsPcaZ;
    double maxAbsDeltaPcaZ;
    double minDaughterPt;
    double minDecayRadius;
    double maxAbsAlpha;
    double maxAbsPairDCA;
    double minDIRA;
    double maxQuality;
    int minNPoints;
  };

  struct HistSet
  {
    TH2F* h_mass_vs_v0pt = nullptr;
    TH3F* h_mass_vs_v0pt_vs_daughterPtMin = nullptr;
    TH3F* h_mass_vs_v0pt_vs_deltaPcaZ = nullptr;
    TH3F* h_mass_vs_v0pt_vs_decayRadius = nullptr;
    TH3F* h_mass_vs_v0pt_vs_alpha = nullptr;
    TH3F* h_mass_vs_v0pt_vs_pairDCA = nullptr;
    TH3F* h_mass_vs_v0pt_vs_DIRA = nullptr;
    TH3F* h_mass_vs_v0pt_vs_npointsMin = nullptr;
    TH3F* h_mass_vs_v0pt_vs_qualityMax = nullptr;
    TH3F* h_mass_vs_v0pt_vs_pcaZ = nullptr;
  };

  bool branchExists(TChain& chain, const char* name)
  {
    return chain.GetBranch(name) != nullptr;
  }

  std::string chargeName(const ChargeCategory category)
  {
    switch (category)
    {
      case ChargeCategory::Unlike: return "unlike";
      case ChargeCategory::Like: return "like";
      case ChargeCategory::PlusPlus: return "plusplus";
      case ChargeCategory::MinusMinus: return "minusminus";
    }
    return "unknown";
  }

  std::vector<ChargeCategory> chargeCategories(double charge1, double charge2)
  {
    std::vector<ChargeCategory> result;
    if (charge1 * charge2 < 0.)
    {
      result.push_back(ChargeCategory::Unlike);
    }
    else if (charge1 * charge2 > 0.)
    {
      result.push_back(ChargeCategory::Like);
      if (charge1 > 0. && charge2 > 0.) result.push_back(ChargeCategory::PlusPlus);
      if (charge1 < 0. && charge2 < 0.) result.push_back(ChargeCategory::MinusMinus);
    }
    return result;
  }

  HistSet bookHistograms(TDirectory* dir,
                         const std::string& selection,
                         const std::string& charge)
  {
    dir->cd();
    const TString tag = TString::Format(" [%s, %s]", selection.c_str(), charge.c_str());

    HistSet h;
    h.h_mass_vs_v0pt = new TH2F(
      "h_mass_Kshort_vs_v0_pt",
      "K^{0}_{S} mass vs V0 p_{T}" + tag + ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]",
      100, 0., 5., 200, 0., 2.);

    h.h_mass_vs_v0pt_vs_daughterPtMin = new TH3F(
      "h_mass_vs_v0pt_vs_daughter_pt_min",
      "Mass vs V0 p_{T} vs lower daughter p_{T}" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];min(p_{T,1},p_{T,2}) [GeV/c]",
      50, 0., 5., 200, 0., 2., 40, 0., 2.);

    h.h_mass_vs_v0pt_vs_deltaPcaZ = new TH3F(
      "h_mass_vs_v0pt_vs_abs_delta_pca_z",
      "Mass vs V0 p_{T} vs |PCA_{z,1}-PCA_{z,2}|" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];|#Delta PCA_{z}| [cm]",
      50, 0., 5., 200, 0., 2., 40, 0., 4.);

    h.h_mass_vs_v0pt_vs_decayRadius = new TH3F(
      "h_mass_vs_v0pt_vs_decay_radius",
      "Mass vs V0 p_{T} vs decay radius" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];decay radius [cm]",
      50, 0., 5., 200, 0., 2., 50, 0., 50.);

    h.h_mass_vs_v0pt_vs_alpha = new TH3F(
      "h_mass_vs_v0pt_vs_abs_alpha",
      "Mass vs V0 p_{T} vs |#alpha|" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];|#alpha|",
      50, 0., 5., 200, 0., 2., 40, 0., 1.);

    h.h_mass_vs_v0pt_vs_pairDCA = new TH3F(
      "h_mass_vs_v0pt_vs_abs_pairDCA",
      "Mass vs V0 p_{T} vs |pair DCA|" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];|pair DCA| [cm]",
      50, 0., 5., 200, 0., 2., 50, 0., 10.);

    h.h_mass_vs_v0pt_vs_DIRA = new TH3F(
      "h_mass_vs_v0pt_vs_DIRA",
      "Mass vs V0 p_{T} vs DIRA" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];DIRA",
      50, 0., 5., 200, 0., 2., 50, 0.5, 1.0);

    h.h_mass_vs_v0pt_vs_npointsMin = new TH3F(
      "h_mass_vs_v0pt_vs_npoints_min",
      "Mass vs V0 p_{T} vs lower daughter N points" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];min(N_{1},N_{2})",
      50, 0., 5., 200, 0., 2., 55, 0., 55.);

    h.h_mass_vs_v0pt_vs_qualityMax = new TH3F(
      "h_mass_vs_v0pt_vs_quality_max",
      "Mass vs V0 p_{T} vs worse daughter quality" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];max(quality_{1},quality_{2})",
      50, 0., 5., 200, 0., 2., 50, 0., 5.);

    h.h_mass_vs_v0pt_vs_pcaZ = new TH3F(
      "h_mass_vs_v0pt_vs_pca_z",
      "Mass vs V0 p_{T} vs candidate PCA z" + tag +
      ";p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}];PCA_{z} [cm]",
      50, 0., 5., 200, 0., 2., 60, -30., 30.);

    return h;
  }

  bool passesSelection(const Selection& s,
                       double pca_z,
                       double absDeltaPcaZ,
                       double daughterPtMin,
                       double decayRadius,
                       double absAlpha,
                       double absPairDCA,
                       double dira,
                       double qualityMax,
                       int npointsMin)
  {
    return std::abs(pca_z) < s.maxAbsPcaZ &&
           absDeltaPcaZ < s.maxAbsDeltaPcaZ &&
           daughterPtMin > s.minDaughterPt &&
           decayRadius > s.minDecayRadius &&
           absAlpha < s.maxAbsAlpha &&
           absPairDCA < s.maxAbsPairDCA &&
           dira > s.minDIRA &&
           qualityMax < s.maxQuality &&
           npointsMin > s.minNPoints;
  }
}

void MakeK0sPairHistograms(const char* inputDir = ".",
                           const char* filePattern = "*.root",
                           const char* outputDir = "output",
                           const char* outputName = "k0s_pair_histograms.root",
                           const char* treeName = "pairTree")
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

  const std::vector<Selection> selections = {
    {"loose",    20.0, 1.50, 0.20, 4.0, 0.90, 6.0, 0.70, 2.0, 20},
    {"baseline", 15.0, 0.50, 0.30, 6.0, 0.80, 4.0, 0.80, 1.0, 30},
    {"tight",    10.0, 0.30, 0.40, 8.0, 0.70, 1.0, 0.90, 0.7, 35}
  };

  gSystem->mkdir(outputDir, kTRUE);
  const TString outputPath = TString::Format("%s/%s", outputDir, outputName);
  std::unique_ptr<TFile> output(TFile::Open(outputPath, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: cannot create " << outputPath << std::endl;
    return;
  }

  std::map<std::string, std::map<std::string, HistSet>> histograms;
  for (const auto& selection : selections)
  {
    TDirectory* selectionDir = output->mkdir(selection.name.c_str());
    for (const auto category : {ChargeCategory::Unlike, ChargeCategory::Like,
                                ChargeCategory::PlusPlus, ChargeCategory::MinusMinus})
    {
      const std::string categoryName = chargeName(category);
      TDirectory* chargeDir = selectionDir->mkdir(categoryName.c_str());
      histograms[selection.name][categoryName] =
        bookHistograms(chargeDir, selection.name, categoryName);
    }
  }

  // Top-level comparison histograms using the user's cuts exactly.
  // Only two charge classes are made here:
  //   unlike sign: charge1*charge2 < 0  (+-)
  //   like sign:   charge1*charge2 > 0  (++ and -- combined)
  output->cd();
  auto makeExactHist = [](const char* name, const char* title) {
    return new TH2F(name, title, 200, 0., 5., 200, 0., 2.);
  };

  TH2F* h_exact1_unlike = makeExactHist(
    "h_exact_cut1_unlike",
    "Exact user cut 1, unlike sign;p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]");
  TH2F* h_exact1_like = makeExactHist(
    "h_exact_cut1_likeSign",
    "Exact user cut 1, like sign (++ + --);p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]");

  TH2F* h_exact2_unlike = makeExactHist(
    "h_exact_cut2_unlike",
    "Exact user cut 2, unlike sign;p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]");
  TH2F* h_exact2_like = makeExactHist(
    "h_exact_cut2_likeSign",
    "Exact user cut 2, like sign (++ + --);p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]");

  TH2F* h_exact3_unlike = makeExactHist(
    "h_exact_cut3_unlike",
    "Exact user cut 3, unlike sign;p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]");
  TH2F* h_exact3_like = makeExactHist(
    "h_exact_cut3_likeSign",
    "Exact user cut 3, like sign (++ + --);p_{T}^{V0} [GeV/c];m_{#pi#pi} [GeV/c^{2}]");

  const Long64_t nEntries = chain.GetEntries();
  std::cout << "Added " << nFiles << " files, entries = " << nEntries << std::endl;

  for (Long64_t entry = 0; entry < nEntries; ++entry)
  {
    chain.GetEntry(entry);
    if (entry % 100000 == 0)
      std::cout << "Processing " << entry << " / " << nEntries << std::endl;

    const auto categories = chargeCategories(charge1, charge2);
    if (categories.empty()) continue;

    const double pt1 = std::hypot(px1, py1);
    const double pt2 = std::hypot(px2, py2);
    const double daughterPtMin = std::min(pt1, pt2);
    const double absDeltaPcaZ = std::abs(pca1_z - pca2_z);
    const double decayRadius = std::hypot(pca_x, pca_y);
    const double absAlpha = std::abs(alpha);
    const double absPairDCA = std::abs(pairDCA);
    const double qualityMax = std::max(quality1, quality2);
    const int npointsMin = std::min(npoints1, npoints2);

    const double pMag = std::sqrt(v0_px*v0_px + v0_py*v0_py + v0_pz*v0_pz);
    const double rMag = std::sqrt(pca_x*pca_x + pca_y*pca_y + pca_z*pca_z);
    const double dira = (pMag > 0. && rMag > 0.)
      ? (v0_px*pca_x + v0_py*pca_y + v0_pz*pca_z)/(pMag*rMag)
      : -2.;

    // Exact user cut 1:
    // pca_z>-15 && pca_z<15 && |pca1_z-pca2_z|<0.5
    // daughter pT>0.3, decay radius>6, |alpha|<0.8, |pairDCA|<1,
    // DIRA>0.8, npoints1/2>30, quality1/2<1.
    const bool exactCut1 =
      pca_z > -15.f && pca_z < 15.f &&
      absDeltaPcaZ < 0.5 &&
      pt1 > 0.3 && pt2 > 0.3 &&
      decayRadius > 6.0 &&
      absAlpha < 0.8 &&
      absPairDCA < 1.0 &&
      dira > 0.8 &&
      npoints1 > 30 && npoints2 > 30 &&
      quality1 < 1.0 && quality2 < 1.0;

    // Exact user cut 2 (the repeated "peak above" cut is identical):
    // -15<pca_z<0, |Delta pca_z|<0.5, daughter pT>0.3,
    // radius>6, |alpha|<0.8, |pairDCA|<4, pz1<0, pz2<0,
    // DIRA>0.8, npoints1/2>30.
    const bool exactCut2 =
      pca_z > -15.f && pca_z < 0.f &&
      absDeltaPcaZ < 0.5 &&
      pt1 > 0.3 && pt2 > 0.3 &&
      decayRadius > 6.0 &&
      absAlpha < 0.8 &&
      absPairDCA < 4.0 &&
      pz1 < 0.f && pz2 < 0.f &&
      dira > 0.8 &&
      npoints1 > 30 && npoints2 > 30;

    // Exact user cut 3: same negative-z cut, but |Delta pca_z|<1.0
    // and no npoints or quality requirement.
    const bool exactCut3 =
      pca_z > -15.f && pca_z < 0.f &&
      absDeltaPcaZ < 1.0 &&
      pt1 > 0.3 && pt2 > 0.3 &&
      decayRadius > 6.0 &&
      absAlpha < 0.8 &&
      absPairDCA < 4.0 &&
      pz1 < 0.f && pz2 < 0.f &&
      dira > 0.8;

    const bool unlikeSign = charge1 * charge2 < 0.f;
    const bool likeSign = charge1 * charge2 > 0.f;

    if (exactCut1)
    {
      if (unlikeSign) h_exact1_unlike->Fill(v0_pt, mass_Kshort);
      if (likeSign) h_exact1_like->Fill(v0_pt, mass_Kshort);
    }
    if (exactCut2)
    {
      if (unlikeSign) h_exact2_unlike->Fill(v0_pt, mass_Kshort);
      if (likeSign) h_exact2_like->Fill(v0_pt, mass_Kshort);
    }
    if (exactCut3)
    {
      if (unlikeSign) h_exact3_unlike->Fill(v0_pt, mass_Kshort);
      if (likeSign) h_exact3_like->Fill(v0_pt, mass_Kshort);
    }

    for (const auto& selection : selections)
    {
      if (!passesSelection(selection, pca_z, absDeltaPcaZ, daughterPtMin,
                           decayRadius, absAlpha, absPairDCA, dira,
                           qualityMax, npointsMin)) continue;

      for (const auto category : categories)
      {
        HistSet& h = histograms.at(selection.name).at(chargeName(category));
        h.h_mass_vs_v0pt->Fill(v0_pt, mass_Kshort);
        h.h_mass_vs_v0pt_vs_daughterPtMin->Fill(v0_pt, mass_Kshort, daughterPtMin);
        h.h_mass_vs_v0pt_vs_deltaPcaZ->Fill(v0_pt, mass_Kshort, absDeltaPcaZ);
        h.h_mass_vs_v0pt_vs_decayRadius->Fill(v0_pt, mass_Kshort, decayRadius);
        h.h_mass_vs_v0pt_vs_alpha->Fill(v0_pt, mass_Kshort, absAlpha);
        h.h_mass_vs_v0pt_vs_pairDCA->Fill(v0_pt, mass_Kshort, absPairDCA);
        h.h_mass_vs_v0pt_vs_DIRA->Fill(v0_pt, mass_Kshort, dira);
        h.h_mass_vs_v0pt_vs_npointsMin->Fill(v0_pt, mass_Kshort, npointsMin);
        h.h_mass_vs_v0pt_vs_qualityMax->Fill(v0_pt, mass_Kshort, qualityMax);
        h.h_mass_vs_v0pt_vs_pcaZ->Fill(v0_pt, mass_Kshort, pca_z);
      }
    }
  }

  output->Write();
  output->Close();
  std::cout << "Wrote output: " << outputPath << std::endl;
}
