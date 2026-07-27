// MakeMineV0PairHistograms.C
//
// One-pass V0 QA for TpcDstV0Finder pairTree.
// Produces 10 cumulative cut levels plus exactCut1/2/3, with:
//   * K0S mass vs V0 pT
//   * Lambda mass vs V0 pT
//   * anti-Lambda mass vs V0 pT
//   * Armenteros-Podolanski qT vs alpha
//   * compact TH3F: V0 pT vs mass vs |pairDCA|
//
// Lambda convention from TpcDstV0Finder:
//   mass_Lambda:     positive daughter = proton, negative daughter = pion
//   mass_AntiLambda: positive daughter = pion,   negative daughter = antiproton
//
// For the Lambda plots, require the proton-hypothesis daughter to have the
// larger transverse momentum. This is a useful background rejection, but it
// is not a fundamental kinematic requirement and can be disabled.
//
// Example:
// root -l -b -q \
// 'MakeMineV0PairHistograms.C("/path/to/V0/files","HITS*_V0.root","output","v0_qa.root","pairTree",true)'

#include <TBranch.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>
#include <TSystem.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  enum class ChargeCategory
  {
    Unlike,
    Like,
    PlusPlus,
    MinusMinus
  };

  struct Selection
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
    int minnpoints;
  };

  struct HistSet
  {
    TH1F* h_k0s_mass = nullptr;
    TH1F* h_lambda_mass = nullptr;
    TH1F* h_antilambda_mass = nullptr;

    TH2F* h_k0s_mass_vs_v0pt = nullptr;
    TH2F* h_lambda_mass_vs_v0pt = nullptr;
    TH2F* h_antilambda_mass_vs_v0pt = nullptr;
    TH2F* h_armenteros_podolanski = nullptr;
    TH2F* h_pair_dca_vs_delta_pca_z = nullptr;

    // Compact 3D histograms.
    TH3F* h3_k0s = nullptr;
    TH3F* h3_lambda = nullptr;
    TH3F* h3_antilambda = nullptr;

    // Daughter pT correlations versus K0S invariant mass.
    TH3F* h3_k0s_pt1_vs_pt2_vs_mass = nullptr;
    TH3F* h3_k0s_ptplus_vs_ptminus_vs_mass = nullptr;
  };
  constexpr double kTwoPi = 2.0 * TMath::Pi();

  double wrapPhi(double phi)
  {
    while (phi >= TMath::Pi()) phi -= kTwoPi;
    while (phi < -TMath::Pi()) phi += kTwoPi;
    return phi;
  }

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

  std::vector<ChargeCategory> chargeCategories(const double charge1,
                                               const double charge2)
  {
    std::vector<ChargeCategory> result;

    if (charge1 * charge2 < 0.)
    {
      result.push_back(ChargeCategory::Unlike);
    }
    else if (charge1 * charge2 > 0.)
    {
      result.push_back(ChargeCategory::Like);

      if (charge1 > 0. && charge2 > 0.)
      {
        result.push_back(ChargeCategory::PlusPlus);
      }
      else if (charge1 < 0. && charge2 < 0.)
      {
        result.push_back(ChargeCategory::MinusMinus);
      }
    }

    return result;
  }

  HistSet bookHistograms(TDirectory* directory,
                         const Selection& selection,
                         const std::string& charge)
  {
    directory->cd();

    const TString tag = TString::Format(
      " [%s, %s]", selection.name.c_str(), charge.c_str());

    HistSet h;

    h.h_k0s_mass = new TH1F(
      "h_mass_Kshort",
      "K^{0}_{S} invariant mass" + tag +
        ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];pairs",
      500, 0.0, 1.0);

    h.h_lambda_mass = new TH1F(
      "h_mass_Lambda",
      "#Lambda invariant mass, p_{T}(p)>p_{T}(#pi)" + tag +
        ";m_{p#pi^{-}} [GeV/c^{2}];pairs",
      300, 1.0, 1.3);

    h.h_antilambda_mass = new TH1F(
      "h_mass_AntiLambda",
      "#bar{#Lambda} invariant mass, p_{T}(#bar{p})>p_{T}(#pi)" + tag +
        ";m_{#bar{p}#pi^{+}} [GeV/c^{2}];pairs",
      300, 1.0, 1.3);

    h.h_k0s_mass_vs_v0pt = new TH2F(
      "h_mass_Kshort_vs_v0_pt",
      "K^{0}_{S} mass vs V0 p_{T}" + tag +
        ";p_{T}^{V0} [GeV/c];m_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
      50, 0., 5.,
      500, 0.0, 1.0);

    h.h_lambda_mass_vs_v0pt = new TH2F(
      "h_mass_Lambda_vs_v0_pt",
      "#Lambda mass vs V0 p_{T}, p_{T}(p)>p_{T}(#pi)" + tag +
        ";p_{T}^{V0} [GeV/c];m_{p#pi^{-}} [GeV/c^{2}]",
      50, 0., 5.,
      300, 1.0, 1.3);

    h.h_antilambda_mass_vs_v0pt = new TH2F(
      "h_mass_AntiLambda_vs_v0_pt",
      "#bar{#Lambda} mass vs V0 p_{T}, p_{T}(#bar{p})>p_{T}(#pi)" + tag +
        ";p_{T}^{V0} [GeV/c];m_{#bar{p}#pi^{+}} [GeV/c^{2}]",
      50, 0., 5.,
      300, 1.0, 1.3);

    h.h_armenteros_podolanski = new TH2F(
      "h_armenteros_podolanski",
      "Armenteros-Podolanski" + tag +
        ";#alpha;q_{T} [GeV/c]",
      240, -1.2, 1.2,
      240, 0.0, 0.40);

    h.h_pair_dca_vs_delta_pca_z = new TH2F(
      "h_pair_dca_vs_abs_delta_pca_z",
      "Pair DCA vs |PCA_{z,1}-PCA_{z,2}|" + tag +
        ";|PCA_{z,1}-PCA_{z,2}| [cm];|pair DCA| [cm]",
      240, 0.0, 2.4,
      240, 0.0, 2.4);

    // Compact 3D binning to control output size:
    //   30 bins in V0 pT, 80 bins in mass, 20 bins in |pairDCA|.
    h.h3_k0s = new TH3F(
      "h3_mass_Kshort_vs_v0_pt_vs_pairDCA",
      "K^{0}_{S}: p_{T} vs mass vs |pair DCA|" + tag +
        ";p_{T}^{V0} [GeV/c];m_{#pi^{+}#pi^{-}} [GeV/c^{2}];|pair DCA| [cm]",
      30, 0.0, 5.0,
      80, 0.30, 0.70,
      20, 0.0, 4.0);

    h.h3_k0s_pt1_vs_pt2_vs_mass = new TH3F(
      "h3_Kshort_pt1_vs_pt2_vs_mass",
      "K^{0}_{S}: daughter p_{T,1} vs p_{T,2} vs mass" + tag +
        ";p_{T,1}^{#pi} [GeV/c];p_{T,2}^{#pi} [GeV/c];"
        "m_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
      50, 0.0, 5.0,
      50, 0.0, 5.0,
      40, 0.40, 0.60);

    h.h3_k0s_ptplus_vs_ptminus_vs_mass = new TH3F(
      "h3_Kshort_ptplus_vs_ptminus_vs_mass",
      "K^{0}_{S}: p_{T}^{#pi^{+}} vs p_{T}^{#pi^{-}} vs mass" + tag +
        ";p_{T}^{#pi^{+}} [GeV/c];p_{T}^{#pi^{-}} [GeV/c];"
        "m_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
      50, 0.0, 5.0,
      50, 0.0, 5.0,
      40, 0.40, 0.60);

    h.h3_lambda = new TH3F(
      "h3_mass_Lambda_vs_v0_pt_vs_pairDCA",
      "#Lambda: p_{T} vs mass vs |pair DCA|" + tag +
        ";p_{T}^{V0} [GeV/c];m_{p#pi^{-}} [GeV/c^{2}];|pair DCA| [cm]",
      30, 0.0, 5.0,
      80, 1.05, 1.25,
      20, 0.0, 4.0);

    h.h3_antilambda = new TH3F(
      "h3_mass_AntiLambda_vs_v0_pt_vs_pairDCA",
      "#bar{#Lambda}: p_{T} vs mass vs |pair DCA|" + tag +
        ";p_{T}^{V0} [GeV/c];m_{#bar{p}#pi^{+}} [GeV/c^{2}];|pair DCA| [cm]",
      30, 0.0, 5.0,
      80, 1.05, 1.25,
      20, 0.0, 4.0);

    return h;
  }

  bool passesSelection(const Selection& selection,
                       const double pcaZ,
                       const double absDeltaPcaZ,
                       const double daughterPtMin,
                       const double decayRadius,
                       const double absAlpha,
                       const double absPairDCA,
                       const double dira,
                       const double qualityMax,
                       const int npointsMin)
  {
    return
      std::abs(pcaZ) < selection.maxAbsPcaZ &&
      absDeltaPcaZ < selection.maxAbsDeltaPcaZ &&
      daughterPtMin > selection.minDaughterPt &&
      decayRadius > selection.minDecayRadius &&
      absAlpha < selection.maxAbsAlpha &&
      absPairDCA < selection.maxAbsPairDCA &&
      dira > selection.minDIRA &&
      qualityMax < selection.maxQuality &&
      npointsMin > selection.minnpoints;
  }
}

void MakeK0sPairHistograms(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputDir = "output",
  const char* outputName = "v0_pair_histograms.root",
  const char* treeName = "pairTree",
  const bool requireProtonHigherPt = true,
  const bool ScaleQualityBy10 = false,
  const Long64_t maxEntries = -1)
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

  const std::vector<std::string> requiredBranches = {
    "mass_Kshort",
    "mass_Lambda",
    "mass_AntiLambda",
    "v0_pt",
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
    "v0_px",
    "v0_py",
    "v0_pz",
    "alpha",
    "qT",
    "pairDCA",
    "npoints1",
    "npoints2",
    "charge1",
    "charge2",
    "quality1",
    "quality2"
  };

  bool missingBranch = false;

  for (const auto& name : requiredBranches)
  {
    if (!branchExists(chain, name.c_str()))
    {
      std::cerr << "ERROR: missing branch "
                << name << std::endl;
      missingBranch = true;
    }
  }

  if (missingBranch)
  {
    return;
  }

  Float_t mass_Kshort = 0.f;
  Float_t mass_Lambda = 0.f;
  Float_t mass_AntiLambda = 0.f;
  Float_t v0_pt = 0.f;

  Float_t pca_x = 0.f;
  Float_t pca_y = 0.f;
  Float_t pca_z = 0.f;
  Float_t pca1_z = 0.f;
  Float_t pca2_z = 0.f;

  Float_t px1 = 0.f;
  Float_t py1 = 0.f;
  Float_t pz1 = 0.f;
  Float_t px2 = 0.f;
  Float_t py2 = 0.f;
  Float_t pz2 = 0.f;

  Float_t v0_px = 0.f;
  Float_t v0_py = 0.f;
  Float_t v0_pz = 0.f;

  Float_t alpha = 0.f;
  Float_t qT = 0.f;
  Float_t pairDCA = 0.f;
  Float_t quality1 = 0.f;
  Float_t quality2 = 0.f;

  Float_t charge1 = 0;
  Float_t charge2 = 0;
  Short_t npoints1 = 0;
  Short_t npoints2 = 0;


  chain.SetBranchAddress("mass_Kshort", &mass_Kshort);
  chain.SetBranchAddress("mass_Lambda", &mass_Lambda);
  chain.SetBranchAddress("mass_AntiLambda", &mass_AntiLambda);
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
  chain.SetBranchAddress("qT", &qT);
  chain.SetBranchAddress("pairDCA", &pairDCA);

  chain.SetBranchAddress("charge1", &charge1);
  chain.SetBranchAddress("charge2", &charge2);
  chain.SetBranchAddress("quality1", &quality1);
  chain.SetBranchAddress("quality2", &quality2);
  chain.SetBranchAddress("npoints1", &npoints1);
  chain.SetBranchAddress("npoints2", &npoints2);

  if (ScaleQualityBy10)
  {
    quality1 *= 10;
    quality2 *= 10;
  }

  // Ten cumulative cut levels.
  //
  // Distances are in cm:
  //   0.50 cm = 5 mm
  //   0.30 cm = 3 mm
  //   0.20 cm = 2 mm
  //   0.10 cm = 1 mm
  //
  // The last levels are intentionally aggressive. Compare signal efficiency
  // and background before choosing a final production selection.
  const std::vector<Selection> selections = {
    {
      "cut00_very_loose",
      "|pca_z|<20, pt>0.20, |pca1_z-pca2_z|<1.50, pairDCA<6.0, DIRA>0.7, Lproj > 2",
      20.0, 1.50, 0.20, 2.0, 0.99, 6.00, 0.70, 20.0, 20
    },
    {
      "cut01_loose",
      "|pca_z|<18, pt>0.20, |pca1_z-pca2_z|<1.00, pairDCA<4.0, DIRA>0.75, Lproj > 2",
      18.0, 1.00, 0.20, 2.0, 0.99, 4.00, 0.75, 20.0, 20
    },
    {
      "cut02_preselection",
      "|pca_z|<15, pt>0.20, |pca1_z-pca2_z|<0.70, pairDCA<3.0, DIRA>0.80, Lproj > 2 cm, n_tpc_clusters>25",
      15.0, 0.70, 0.20, 2.0, 0.99, 3.00, 0.80, 18.0, 25
    },
    {
      "cut03_baseline",
      "|pca_z|<15, pt>0.20, |pca1_z-pca2_z|<0.50, pairDCA<2.0, DIRA>0.85, Lproj > 2 cm, n_tpc_hits>30",
      15.0, 0.50, 0.20, 2.0, 0.99, 2.00, 0.85, 15.0, 30
    },
    {
      "cut04_pairDCA_15mm",
      "|z|<12, pt>0.30, dz<0.50, pairDCA<1.5",
      12.0, 0.50, 0.20, 2.0, 0.99, 1.50, 0.88, 15.0, 30
    },
    {
      "cut05_pairDCA_10mm",
      "|z|<10, pt>0.30, dz<0.30, pairDCA<1.0",
      10.0, 0.40, 0.20, 2.0, 0.99, 1.00, 0.90, 13.0, 30
    },
    {
      "cut06_pairDCA_7mm",
      "|z|<10, pt>0.3, dz<0.25, pairDCA<0.7",
      10.0, 0.25, 0.20, 2.0, 0.99, 0.70, 0.93, 12.0, 32
    },
    {
      "cut07_pairDCA_5mm",
      "|z|<10, pt>0.30, dz<0.20, pairDCA<0.5",
      10.0, 0.20, 0.20, 2.0, 0.99, 0.50, 0.95, 10.0, 32
    },
    {
      "cut08_pairDCA_3mm",
      "|z|<8, pt>0.30, dz<0.15, pairDCA<0.3",
      8.00, 0.15, 0.20, 2.0, 0.99, 0.30, 0.98, 8.0, 35
    },
    {
      "cut09_pairDCA_2mm",
      "|z|<8, pt>0.20, dz<0.10, pairDCA<0.2",
      8.00, 0.10, 0.20, 2.0, 0.99, 0.20, 0.99, 6.0, 35
    }
  };

  gSystem->mkdir(outputDir, kTRUE);

  const TString outputPath =
    TString::Format("%s/%s", outputDir, outputName);

  std::unique_ptr<TFile> output(
    TFile::Open(outputPath, "RECREATE"));

  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: cannot create "
              << outputPath << std::endl;
    return;
  }

  output->cd();

  TH1I* h_cutflow_unlike = new TH1I(
    "h_cutflow_unlike",
    "Unlike-sign cumulative cut flow;selection;pair count",
    selections.size(), 0, selections.size());

  TH1I* h_cutflow_like = new TH1I(
    "h_cutflow_like",
    "Like-sign cumulative cut flow;selection;pair count",
    selections.size(), 0, selections.size());

  for (std::size_t i = 0; i < selections.size(); ++i)
  {
    h_cutflow_unlike->GetXaxis()->SetBinLabel(
      i + 1, selections[i].name.c_str());

    h_cutflow_like->GetXaxis()->SetBinLabel(
      i + 1, selections[i].name.c_str());
  }

  std::map<std::string, std::map<std::string, HistSet>> histograms;

  for (const auto& selection : selections)
  {
    TDirectory* selectionDir =
      output->mkdir(selection.name.c_str());

    selectionDir->cd();

    TNamed selectionInfo(
      "selection",
      selection.description.c_str());
    selectionInfo.Write();

    for (const auto category : {
           ChargeCategory::Unlike,
           ChargeCategory::Like,
           ChargeCategory::PlusPlus,
           ChargeCategory::MinusMinus})
    {
      const std::string categoryName =
        chargeName(category);

      TDirectory* chargeDir =
        selectionDir->mkdir(categoryName.c_str());

      histograms[selection.name][categoryName] =
        bookHistograms(chargeDir, selection, categoryName);
    }
  }

  // Separate exact user selections. These do not replace the 10-cut scan.
  const std::vector<std::pair<std::string, std::string>> exactSelections = {
    {
      "exactCut1",
      "-15<pca_z<15, |Delta pca_z|<0.5, daughter pT>0.3, "
      "radius>6, |alpha|<0.8, |pairDCA|<1, DIRA>0.8, "
      "npoints>30, quality<20"
    },
    {
      "exactCut2",
      "-15<pca_z<0, |Delta pca_z|<0.5, daughter pT>0.3, "
      "radius>6, |alpha|<0.8, |pairDCA|<4, pz1/2<0, "
      "DIRA>0.8, npoints>30"
    },
    {
      "exactCut3",
      "-15<pca_z<0, |Delta pca_z|<1.0, daughter pT>0.3, "
      "radius>6, |alpha|<0.8, |pairDCA|<4, pz1/2<0, DIRA>0.8"
    }
  };

  std::map<std::string, std::map<std::string, HistSet>> exactHistograms;

  for (const auto& exactSelection : exactSelections)
  {
    TDirectory* selectionDir =
      output->mkdir(exactSelection.first.c_str());

    selectionDir->cd();

    TNamed selectionInfo(
      "selection",
      exactSelection.second.c_str());
    selectionInfo.Write();

    Selection bookingSelection{
      exactSelection.first, exactSelection.second,
      0., 0., 0., 0., 0., 0., 0., 0., 0
    };

    for (const auto category : {
           ChargeCategory::Unlike,
           ChargeCategory::Like,
           ChargeCategory::PlusPlus,
           ChargeCategory::MinusMinus})
    {
      const std::string categoryName = chargeName(category);
      TDirectory* chargeDir =
        selectionDir->mkdir(categoryName.c_str());

      exactHistograms[exactSelection.first][categoryName] =
        bookHistograms(chargeDir, bookingSelection, categoryName);
    }
  }

  const Long64_t totalEntries = chain.GetEntries();

  const Long64_t entriesToProcess =
    (maxEntries >= 0)
      ? std::min(totalEntries, maxEntries)
      : totalEntries;

  std::cout
    << "Added " << nFiles
    << " files, total entries = " << totalEntries
    << ", processing = " << entriesToProcess
    << std::endl;

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


    
    const double phiPos = std::atan2(charge1==1 ? py1 : py2, charge1==1 ? px1 : px2);
    const double phiNeg = std::atan2(charge1==1 ? py2 : py1, charge1==1 ? px2 : px1);
    const double deltaPhi = wrapPhi(phiPos - phiNeg);
    if (deltaPhi <0.8 - 0.4 * ( v0_pt  < 2.0 ? v0_pt : 2.0 ) )
    {
      continue;
    }
    const auto categories =
      chargeCategories(charge1, charge2);

    if (categories.empty())
    {
      continue;
    }

    const double pt1 = std::hypot(px1, py1);
    const double pt2 = std::hypot(px2, py2);

    const double daughterPtMin =
      std::min(pt1, pt2);

    const double absDeltaPcaZ =
      std::abs(pca1_z - pca2_z);

    const double decayRadius =
      std::hypot(pca_x, pca_y);

    const double absAlpha =
      std::abs(alpha);

    const double absPairDCA =
      std::abs(pairDCA);

    const double qualityMax =
      std::max(quality1, quality2);

    const int npointsMin =
      std::min<int>(npoints1, npoints2);

    const double v0Momentum =
      std::sqrt(
        v0_px * v0_px +
        v0_py * v0_py +
        v0_pz * v0_pz);

    const double flightLength =
      std::sqrt(
        pca_x * pca_x +
        pca_y * pca_y +
        pca_z * pca_z);

    const double dira =
      (v0Momentum > 0. && flightLength > 0.)
        ? (v0_px * pca_x +
           v0_py * pca_y +
           v0_pz * pca_z) /
          (v0Momentum * flightLength)
        : -2.;

    const bool unlikeSign =
      charge1 * charge2 < 0;

    const bool likeSign =
      charge1 * charge2 > 0;

    double positivePt = -1.;
    double negativePt = -1.;

    if (unlikeSign)
    {
      if (charge1 > 0)
      {
        positivePt = pt1;
        negativePt = pt2;
      }
      else
      {
        positivePt = pt2;
        negativePt = pt1;
      }
    }

    const int effectiveNpoints1 = npoints1;

    const int effectiveNpoints2 = npoints2;

    const bool exactCut1 =
      pca_z > -15.f && pca_z < 15.f &&
      absDeltaPcaZ < 0.5 &&
      pt1 > 0.3 && pt2 > 0.3 &&
      decayRadius > 6.0 &&
      absAlpha < 0.8 &&
      absPairDCA < 1.0 &&
      dira > 0.8 &&
      effectiveNpoints1 > 30 && effectiveNpoints2 > 30 &&
      quality1 < 20.0 && quality2 < 20.0;

    const bool exactCut2 =
      pca_z > -15.f && pca_z < 0.f &&
      absDeltaPcaZ < 0.5 &&
      pt1 > 0.3 && pt2 > 0.3 &&
      decayRadius > 6.0 &&
      absAlpha < 0.8 &&
      absPairDCA < 4.0 &&
      pz1 < 0.f && pz2 < 0.f &&
      dira > 0.8 &&
      effectiveNpoints1 > 30 && effectiveNpoints2 > 30;

    const bool exactCut3 =
      pca_z > -15.f && pca_z < 0.f &&
      absDeltaPcaZ < 1.0 &&
      pt1 > 0.3 && pt2 > 0.3 &&
      decayRadius > 6.0 &&
      absAlpha < 0.8 &&
      absPairDCA < 4.0 &&
      pz1 < 0.f && pz2 < 0.f &&
      dira > 0.8;

    // Lambda: positive daughter carries the proton hypothesis.
    const bool lambdaProtonPtPass =
      unlikeSign &&
      (!requireProtonHigherPt ||
       positivePt > negativePt);

    // anti-Lambda: negative daughter carries the antiproton hypothesis.
    const bool antiLambdaProtonPtPass =
      unlikeSign &&
      (!requireProtonHigherPt ||
       negativePt > positivePt);

    for (std::size_t selectionIndex = 0;
         selectionIndex < selections.size();
         ++selectionIndex)
    {
      const auto& selection =
        selections[selectionIndex];

      if (!passesSelection(
            selection,
            pca_z,
            absDeltaPcaZ,
            daughterPtMin,
            decayRadius,
            absAlpha,
            absPairDCA,
            dira,
            qualityMax,
            npointsMin))
      {
        continue;
      }

      if (unlikeSign)
      {
        h_cutflow_unlike->Fill(selectionIndex + 0.5);
      }

      if (likeSign)
      {
        h_cutflow_like->Fill(selectionIndex + 0.5);
      }

      for (const auto category : categories)
      {
        HistSet& h =
          histograms
            .at(selection.name)
            .at(chargeName(category));

        h.h_k0s_mass->Fill(mass_Kshort);

        h.h_k0s_mass_vs_v0pt->Fill(
          v0_pt, mass_Kshort);

        h.h3_k0s->Fill(
          v0_pt, mass_Kshort, absPairDCA);

        h.h3_k0s_pt1_vs_pt2_vs_mass->Fill(
          pt1, pt2, mass_Kshort);

        if (category == ChargeCategory::Unlike)
        {
          h.h3_k0s_ptplus_vs_ptminus_vs_mass->Fill(
            positivePt, negativePt, mass_Kshort);
        }

        h.h_armenteros_podolanski->Fill(
          alpha, qT);

        h.h_pair_dca_vs_delta_pca_z->Fill(
          absDeltaPcaZ, absPairDCA);

        // Lambda and anti-Lambda are physically meaningful only
        // for unlike-sign pairs.
        if (category == ChargeCategory::Unlike)
        {
          if (lambdaProtonPtPass)
          {
            h.h_lambda_mass->Fill(mass_Lambda);
            h.h_lambda_mass_vs_v0pt->Fill(
              v0_pt, mass_Lambda);

            h.h3_lambda->Fill(
              v0_pt, mass_Lambda, absPairDCA);
          }

          if (antiLambdaProtonPtPass)
          {
            h.h_antilambda_mass->Fill(mass_AntiLambda);
            h.h_antilambda_mass_vs_v0pt->Fill(
              v0_pt, mass_AntiLambda);

            h.h3_antilambda->Fill(
              v0_pt, mass_AntiLambda, absPairDCA);
          }
        }
      }
    }

    const std::vector<std::pair<std::string, bool>> exactPasses = {
      {"exactCut1", exactCut1},
      {"exactCut2", exactCut2},
      {"exactCut3", exactCut3}
    };

    for (const auto& exactPass : exactPasses)
    {
      if (!exactPass.second)
      {
        continue;
      }

      for (const auto category : categories)
      {
        HistSet& h =
          exactHistograms
            .at(exactPass.first)
            .at(chargeName(category));

        h.h_k0s_mass->Fill(mass_Kshort);
        h.h_k0s_mass_vs_v0pt->Fill(v0_pt, mass_Kshort);
        h.h3_k0s->Fill(v0_pt, mass_Kshort, absPairDCA);
        h.h3_k0s_pt1_vs_pt2_vs_mass->Fill(
          pt1, pt2, mass_Kshort);

        if (category == ChargeCategory::Unlike)
        {
          h.h3_k0s_ptplus_vs_ptminus_vs_mass->Fill(
            positivePt, negativePt, mass_Kshort);
        }

        h.h_armenteros_podolanski->Fill(alpha, qT);
        h.h_pair_dca_vs_delta_pca_z->Fill(
          absDeltaPcaZ, absPairDCA);

        if (category == ChargeCategory::Unlike)
        {
          if (lambdaProtonPtPass)
          {
            h.h_lambda_mass->Fill(mass_Lambda);
            h.h_lambda_mass_vs_v0pt->Fill(v0_pt, mass_Lambda);
            h.h3_lambda->Fill(v0_pt, mass_Lambda, absPairDCA);
          }

          if (antiLambdaProtonPtPass)
          {
            h.h_antilambda_mass->Fill(mass_AntiLambda);
            h.h_antilambda_mass_vs_v0pt->Fill(v0_pt, mass_AntiLambda);
            h.h3_antilambda->Fill(
              v0_pt, mass_AntiLambda, absPairDCA);
          }
        }
      }
    }
  }

  output->Write();
  output->Close();

  std::cout
    << "Wrote output: "
    << outputPath
    << std::endl;
}
