// MakeMineV0PairHistograms.C
//
// One-pass V0 QA for TpcDstV0Finder pairTree.
// Produces 10 cumulative cut levels, with:
//   * K0S mass vs V0 pT
//   * Lambda mass vs V0 pT
//   * anti-Lambda mass vs V0 pT
//   * Armenteros-Podolanski qT vs alpha
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
#include <TH2F.h>
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
    int minNclusters;
  };

  struct HistSet
  {
    TH2F* h_k0s_mass_vs_v0pt = nullptr;
    TH2F* h_lambda_mass_vs_v0pt = nullptr;
    TH2F* h_antilambda_mass_vs_v0pt = nullptr;
    TH2F* h_armenteros_podolanski = nullptr;
    TH2F* h_pair_dca_vs_delta_pca_z = nullptr;
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

    h.h_k0s_mass_vs_v0pt = new TH2F(
      "h_mass_Kshort_vs_v0_pt",
      "K^{0}_{S} mass vs V0 p_{T}" + tag +
        ";p_{T}^{V0} [GeV/c];m_{#pi^{+}#pi^{-}} [GeV/c^{2}]",
      100, 0., 5.,
      400, 0.30, 0.70);

    h.h_lambda_mass_vs_v0pt = new TH2F(
      "h_mass_Lambda_vs_v0_pt",
      "#Lambda mass vs V0 p_{T}, p_{T}(p)>p_{T}(#pi)" + tag +
        ";p_{T}^{V0} [GeV/c];m_{p#pi^{-}} [GeV/c^{2}]",
      100, 0., 5.,
      400, 1.05, 1.25);

    h.h_antilambda_mass_vs_v0pt = new TH2F(
      "h_mass_AntiLambda_vs_v0_pt",
      "#bar{#Lambda} mass vs V0 p_{T}, p_{T}(#bar{p})>p_{T}(#pi)" + tag +
        ";p_{T}^{V0} [GeV/c];m_{#bar{p}#pi^{+}} [GeV/c^{2}]",
      100, 0., 5.,
      400, 1.05, 1.25);

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
                       const int nclustersMin)
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
      nclustersMin > selection.minNclusters;
  }
}

void MakeMineV0PairHistograms(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputDir = "output",
  const char* outputName = "v0_pair_histograms.root",
  const char* treeName = "pairTree",
  const bool requireProtonHigherPt = true,
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
    "nclusters1",
    "nclusters2",
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

  Short_t charge1 = 0;
  Short_t charge2 = 0;
  UShort_t nclusters1 = 0;
  UShort_t nclusters2 = 0;

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
  chain.SetBranchAddress("nclusters1", &nclusters1);
  chain.SetBranchAddress("nclusters2", &nclusters2);

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
      "|z|<20, dz<1.50, pairDCA<6.0",
      20.0, 1.50, 0.20, 4.0, 0.95, 6.00, 0.70, 20.0, 20
    },
    {
      "cut01_loose",
      "|z|<18, dz<1.00, pairDCA<4.0",
      18.0, 1.00, 0.20, 4.0, 0.90, 4.00, 0.75, 20.0, 20
    },
    {
      "cut02_preselection",
      "|z|<15, dz<0.70, pairDCA<3.0",
      15.0, 0.70, 0.25, 5.0, 0.85, 3.00, 0.80, 18.0, 25
    },
    {
      "cut03_baseline",
      "|z|<15, dz<0.50, pairDCA<2.0",
      15.0, 0.50, 0.30, 6.0, 0.80, 2.00, 0.85, 15.0, 30
    },
    {
      "cut04_pairDCA_15mm",
      "|z|<12, dz<0.40, pairDCA<1.5",
      12.0, 0.40, 0.30, 6.0, 0.80, 1.50, 0.88, 14.0, 30
    },
    {
      "cut05_pairDCA_10mm",
      "|z|<10, dz<0.30, pairDCA<1.0",
      10.0, 0.30, 0.35, 6.0, 0.75, 1.00, 0.90, 13.0, 30
    },
    {
      "cut06_pairDCA_7mm",
      "|z|<10, dz<0.25, pairDCA<0.7",
      10.0, 0.25, 0.35, 7.0, 0.75, 0.70, 0.93, 12.0, 32
    },
    {
      "cut07_pairDCA_5mm",
      "|z|<10, dz<0.20, pairDCA<0.5",
      10.0, 0.20, 0.40, 7.0, 0.70, 0.50, 0.95, 10.0, 32
    },
    {
      "cut08_pairDCA_3mm",
      "|z|<8, dz<0.15, pairDCA<0.3",
      8.0, 0.15, 0.40, 8.0, 0.70, 0.30, 0.98, 8.0, 35
    },
    {
      "cut09_pairDCA_2mm",
      "|z|<8, dz<0.10, pairDCA<0.2",
      8.0, 0.10, 0.45, 8.0, 0.65, 0.20, 0.995, 6.0, 35
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

    const int nclustersMin =
      std::min<int>(nclusters1, nclusters2);

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
            nclustersMin))
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

        h.h_k0s_mass_vs_v0pt->Fill(
          v0_pt, mass_Kshort);

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
            h.h_lambda_mass_vs_v0pt->Fill(
              v0_pt, mass_Lambda);
          }

          if (antiLambdaProtonPtPass)
          {
            h.h_antilambda_mass_vs_v0pt->Fill(
              v0_pt, mass_AntiLambda);
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
