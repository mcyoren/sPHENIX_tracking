// MakeMineV0PairHistograms.C
//
// One-pass V0 QA for TpcDstV0Finder pairTree.
// Produces 10 cumulative cut levels plus exactCut1/2/3, with:
//   * K0S mass vs V0 pT
//   * Lambda mass vs V0 pT
//   * anti-Lambda mass vs V0 pT
//   * phi -> K+K- mass QA with prompt-track cuts
//   * D0/anti-D0 -> K pi mass QA with prompt-track cuts
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
#include <TNamed.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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

    // Negative bound disables that side of the dE/dx selection.
    double pionDedxMin{-1.0};
    double pionDedxMax{400.0};
    double protonDedxMin{-1.0};
    double protonDedxMax{-1.0};
  };

  struct PromptSelection
  {
    std::string name;
    std::string description;

    double minDaughterPt;
    unsigned int minNtpcClusters;
    double maxQuality;
    double maxAbsTrackDcaXY;
    double maxAbsTrackDcaZ;
    double maxAbsPrimaryPcaZ;
    double maxAbsDeltaPrimaryPcaZ;
    double maxAbsPromptPairDCA;

    // Negative bounds disable the corresponding dE/dx side.
    double pionDedxMin{-1.0};
    double pionDedxMax{-1.0};
    double kaonDedxMin{-1.0};
    double kaonDedxMax{-1.0};
  };

  struct HistSet
  {
    TH1F* h_k0s_mass = nullptr;
    TH1F* h_lambda_mass = nullptr;
    TH1F* h_antilambda_mass = nullptr;
    TH1F* h_phi_mass = nullptr;
    TH1F* h_d0_mass = nullptr;
    TH1F* h_antid0_mass = nullptr;

    TH2F* h_k0s_mass_vs_v0pt = nullptr;
    TH2F* h_lambda_mass_vs_v0pt = nullptr;
    TH2F* h_antilambda_mass_vs_v0pt = nullptr;
    TH2F* h_phi_mass_vs_v0pt = nullptr;
    TH2F* h_d0_mass_vs_v0pt = nullptr;
    TH2F* h_antid0_mass_vs_v0pt = nullptr;
    TH2F* h_armenteros_podolanski = nullptr;
    TH2F* h_pair_dca_vs_delta_pca_z = nullptr;

    // Prompt-primary QA. Filled only in the promptMesons directories.
    TH1F* h_primary_pca_dz = nullptr;
    TH2F* h_primary_pca_z1_vs_z2 = nullptr;
    TH2F* h_max_dca_xy_vs_primary_pca_dz = nullptr;

    // Compact 3D histograms.
    TH3F* h3_k0s = nullptr;
    TH3F* h3_lambda = nullptr;
    TH3F* h3_antilambda = nullptr;
    TH3F* h3_phi = nullptr;
    TH3F* h3_d0 = nullptr;
    TH3F* h3_antid0 = nullptr;

    // Daughter pT correlations versus K0S invariant mass.
    TH3F* h3_k0s_pt1_vs_pt2_vs_mass = nullptr;
    TH3F* h3_k0s_ptplus_vs_ptminus_vs_mass = nullptr;
  };
  constexpr double kTwoPi = 2.0 * TMath::Pi();
  constexpr double kPionMass = 0.13957039;
  constexpr double kKaonMass = 0.493677;

  double invariantMass(const double px1,
                       const double py1,
                       const double pz1,
                       const double mass1,
                       const double px2,
                       const double py2,
                       const double pz2,
                       const double mass2)
  {
    const double p1Squared = px1 * px1 + py1 * py1 + pz1 * pz1;
    const double p2Squared = px2 * px2 + py2 * py2 + pz2 * pz2;

    const double energy1 = std::sqrt(p1Squared + mass1 * mass1);
    const double energy2 = std::sqrt(p2Squared + mass2 * mass2);

    const double totalPx = px1 + px2;
    const double totalPy = py1 + py2;
    const double totalPz = pz1 + pz2;
    const double totalEnergy = energy1 + energy2;

    const double massSquared =
      totalEnergy * totalEnergy -
      totalPx * totalPx -
      totalPy * totalPy -
      totalPz * totalPz;

    return std::sqrt(std::max(0.0, massSquared));
  }

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

  bool passesDedx(const double value,
                  const double minValue,
                  const double maxValue)
  {
    return
      (minValue < 0.0 || value > minValue) &&
      (maxValue < 0.0 || value < maxValue);
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

    h.h_phi_mass = new TH1F(
      "h_mass_Phi",
      "#phi invariant mass" + tag +
        ";m_{K^{+}K^{-}} [GeV/c^{2}];pairs",
      300, 0.95, 1.10);

    h.h_d0_mass = new TH1F(
      "h_mass_D0",
      "D^{0} invariant mass" + tag +
        ";m_{K^{-}#pi^{+}} [GeV/c^{2}];pairs",
      350, 1.65, 2.10);

    h.h_antid0_mass = new TH1F(
      "h_mass_AntiD0",
      "#bar{D}^{0} invariant mass" + tag +
        ";m_{K^{+}#pi^{-}} [GeV/c^{2}];pairs",
      350, 1.65, 2.10);

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

    h.h_phi_mass_vs_v0pt = new TH2F(
      "h_mass_Phi_vs_v0_pt",
      "#phi mass vs pair p_{T}" + tag +
        ";p_{T}^{pair} [GeV/c];m_{K^{+}K^{-}} [GeV/c^{2}]",
      50, 0., 5.,
      300, 0.95, 1.10);

    h.h_d0_mass_vs_v0pt = new TH2F(
      "h_mass_D0_vs_v0_pt",
      "D^{0} mass vs pair p_{T}" + tag +
        ";p_{T}^{pair} [GeV/c];m_{K^{-}#pi^{+}} [GeV/c^{2}]",
      50, 0., 5.,
      350, 1.65, 2.10);

    h.h_antid0_mass_vs_v0pt = new TH2F(
      "h_mass_AntiD0_vs_v0_pt",
      "#bar{D}^{0} mass vs pair p_{T}" + tag +
        ";p_{T}^{pair} [GeV/c];m_{K^{+}#pi^{-}} [GeV/c^{2}]",
      50, 0., 5.,
      350, 1.65, 2.10);

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

    h.h_primary_pca_dz = new TH1F(
      "h_primary_pca_dz",
      "Primary beam-axis PCA z difference" + tag +
        ";|z_{PCA,1}^{primary}-z_{PCA,2}^{primary}| [cm];pairs",
      240, 0.0, 2.4);

    h.h_primary_pca_z1_vs_z2 = new TH2F(
      "h_primary_pca_z1_vs_z2",
      "Primary beam-axis PCA z correlation" + tag +
        ";z_{PCA,1}^{primary} [cm];z_{PCA,2}^{primary} [cm]",
      160, -40.0, 40.0,
      160, -40.0, 40.0);

    h.h_max_dca_xy_vs_primary_pca_dz = new TH2F(
      "h_max_dca_xy_vs_primary_pca_dz",
      "Track DCA_{xy} vs primary PCA z difference" + tag +
        ";|z_{PCA,1}^{primary}-z_{PCA,2}^{primary}| [cm];max |DCA_{xy}| [cm]",
      240, 0.0, 2.4,
      200, 0.0, 5.0);

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

    h.h3_phi = new TH3F(
      "h3_mass_Phi_vs_v0_pt_vs_pairDCA",
      "#phi: p_{T} vs mass vs |pair DCA|" + tag +
        ";p_{T}^{pair} [GeV/c];m_{K^{+}K^{-}} [GeV/c^{2}];|pair DCA| [cm]",
      30, 0.0, 5.0,
      80, 0.98, 1.06,
      20, 0.0, 2.0);

    h.h3_d0 = new TH3F(
      "h3_mass_D0_vs_v0_pt_vs_pairDCA",
      "D^{0}: p_{T} vs mass vs |pair DCA|" + tag +
        ";p_{T}^{pair} [GeV/c];m_{K^{-}#pi^{+}} [GeV/c^{2}];|pair DCA| [cm]",
      30, 0.0, 5.0,
      90, 1.70, 2.05,
      20, 0.0, 2.0);

    h.h3_antid0 = new TH3F(
      "h3_mass_AntiD0_vs_v0_pt_vs_pairDCA",
      "#bar{D}^{0}: p_{T} vs mass vs |pair DCA|" + tag +
        ";p_{T}^{pair} [GeV/c];m_{K^{+}#pi^{-}} [GeV/c^{2}];|pair DCA| [cm]",
      30, 0.0, 5.0,
      90, 1.70, 2.05,
      20, 0.0, 2.0);

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


  bool passesPromptSelection(const PromptSelection& selection,
                             const double daughterPtMin,
                             const unsigned int ntpcClustersMin,
                             const double qualityMax,
                             const double maxAbsTrackDcaXY,
                             const double maxAbsTrackDcaZ,
                             const double maxAbsPrimaryPcaZ,
                             const double absDeltaPrimaryPcaZ,
                             const double absPromptPairDCA)
  {
    return
      daughterPtMin > selection.minDaughterPt &&
      ntpcClustersMin >= selection.minNtpcClusters &&
      (selection.maxQuality < 0.0 ||
       qualityMax < selection.maxQuality) &&
      (selection.maxAbsTrackDcaXY < 0.0 ||
       maxAbsTrackDcaXY < selection.maxAbsTrackDcaXY) &&
      (selection.maxAbsTrackDcaZ < 0.0 ||
       maxAbsTrackDcaZ < selection.maxAbsTrackDcaZ) &&
      (selection.maxAbsPrimaryPcaZ < 0.0 ||
       maxAbsPrimaryPcaZ < selection.maxAbsPrimaryPcaZ) &&
      (selection.maxAbsDeltaPrimaryPcaZ < 0.0 ||
       absDeltaPrimaryPcaZ < selection.maxAbsDeltaPrimaryPcaZ) &&
      (selection.maxAbsPromptPairDCA < 0.0 ||
       (std::isfinite(absPromptPairDCA) &&
        absPromptPairDCA < selection.maxAbsPromptPairDCA));
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
  const Long64_t maxEntries = -1,
  const bool usePrimaryVertexKinematicsForPrompt = true,
  const double beamX = 0.158,
  const double beamY = 0.285,
  const double beamZ = 0.0)
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
    "Lproj",
    "dca_xy1",
    "dca_z1",
    "dca_xy2",
    "dca_z2",
    "npoints1",
    "npoints2",
    "charge1",
    "charge2",
    "quality1",
    "quality2",
    "dedx_1",
    "dedx_2"
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

  if (usePrimaryVertexKinematicsForPrompt)
  {
    const std::vector<std::string> requiredPrimaryBranches = {
      "primary_px1",
      "primary_py1",
      "primary_pz1",
      "primary_px2",
      "primary_py2",
      "primary_pz2",
      "primary_pca1_z",
      "primary_pca2_z",
      "primary_pca_dz",
      "primary_pca_valid",
      "prompt_pairDCA",
      "prompt_pca_valid",
      "ntpc_clusters1",
      "ntpc_clusters2"
    };

    for (const auto& name : requiredPrimaryBranches)
    {
      if (!branchExists(chain, name.c_str()))
      {
        std::cerr
          << "ERROR: usePrimaryVertexKinematicsForPrompt=true, "
          << "but branch " << name << " is missing" << std::endl;
        missingBranch = true;
      }
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

  Float_t primary_px1 = 0.f;
  Float_t primary_py1 = 0.f;
  Float_t primary_pz1 = 0.f;
  Float_t primary_px2 = 0.f;
  Float_t primary_py2 = 0.f;
  Float_t primary_pz2 = 0.f;
  Float_t primary_pca1_z = 0.f;
  Float_t primary_pca2_z = 0.f;
  Float_t primary_pca_dz = 0.f;
  Int_t primary_pca_valid = 0;
  Float_t prompt_pairDCA = 0.f;
  Int_t prompt_pca_valid = 0;

  Float_t v0_px = 0.f;
  Float_t v0_py = 0.f;
  Float_t v0_pz = 0.f;

  Float_t alpha = 0.f;
  Float_t qT = 0.f;
  Float_t pairDCA = 0.f;
  Float_t Lproj = 0.f;
  Float_t dca_xy1 = 0.f;
  Float_t dca_z1 = 0.f;
  Float_t dca_xy2 = 0.f;
  Float_t dca_z2 = 0.f;
  Float_t quality1 = 0.f;
  Float_t quality2 = 0.f;
  Float_t dedx_1 = 0.f;
  Float_t dedx_2 = 0.f;

  Float_t charge1 = 0;
  Float_t charge2 = 0;
  Short_t npoints1 = 0;
  Short_t npoints2 = 0;
  UInt_t ntpc_clusters1 = 0;
  UInt_t ntpc_clusters2 = 0;


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

  if (usePrimaryVertexKinematicsForPrompt)
  {
    chain.SetBranchAddress("primary_px1", &primary_px1);
    chain.SetBranchAddress("primary_py1", &primary_py1);
    chain.SetBranchAddress("primary_pz1", &primary_pz1);
    chain.SetBranchAddress("primary_px2", &primary_px2);
    chain.SetBranchAddress("primary_py2", &primary_py2);
    chain.SetBranchAddress("primary_pz2", &primary_pz2);
    chain.SetBranchAddress("primary_pca1_z", &primary_pca1_z);
    chain.SetBranchAddress("primary_pca2_z", &primary_pca2_z);
    chain.SetBranchAddress("primary_pca_dz", &primary_pca_dz);
    chain.SetBranchAddress("primary_pca_valid", &primary_pca_valid);
    chain.SetBranchAddress("prompt_pairDCA", &prompt_pairDCA);
    chain.SetBranchAddress("prompt_pca_valid", &prompt_pca_valid);
    chain.SetBranchAddress("ntpc_clusters1", &ntpc_clusters1);
    chain.SetBranchAddress("ntpc_clusters2", &ntpc_clusters2);
  }

  chain.SetBranchAddress("v0_px", &v0_px);
  chain.SetBranchAddress("v0_py", &v0_py);
  chain.SetBranchAddress("v0_pz", &v0_pz);

  chain.SetBranchAddress("alpha", &alpha);
  chain.SetBranchAddress("qT", &qT);
  chain.SetBranchAddress("pairDCA", &pairDCA);
  chain.SetBranchAddress("Lproj", &Lproj);
  chain.SetBranchAddress("dca_xy1", &dca_xy1);
  chain.SetBranchAddress("dca_z1", &dca_z1);
  chain.SetBranchAddress("dca_xy2", &dca_xy2);
  chain.SetBranchAddress("dca_z2", &dca_z2);

  chain.SetBranchAddress("charge1", &charge1);
  chain.SetBranchAddress("charge2", &charge2);
  chain.SetBranchAddress("quality1", &quality1);
  chain.SetBranchAddress("quality2", &quality2);
  chain.SetBranchAddress("dedx_1", &dedx_1);
  chain.SetBranchAddress("dedx_2", &dedx_2);
  chain.SetBranchAddress("npoints1", &npoints1);
  chain.SetBranchAddress("npoints2", &npoints2);

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
      12.0, 0.50, 0.30, 2.0, 0.99, 1.50, 0.88, 15.0, 30
    },
    {
      "cut05_pairDCA_10mm",
      "|z|<10, pt>0.30, dz<0.30, pairDCA<1.0",
      10.0, 0.30, 0.30, 2.0, 0.99, 1.00, 0.90, 13.0, 30
    },
    {
      "cut06_pairDCA_7mm",
      "|z|<10, pt>0.3, dz<0.25, pairDCA<0.7",
      10.0, 0.25, 0.30, 2.0, 0.99, 0.70, 0.93, 12.0, 32
    },
    {
      "cut07_pairDCA_5mm",
      "|z|<10, pt>0.30, dz<0.20, pairDCA<0.5",
      10.0, 0.20, 0.30, 2.0, 0.99, 0.50, 0.95, 10.0, 32
    },
    {
      "cut08_pairDCA_3mm",
      "|z|<8, pt>0.30, dz<0.15, pairDCA<0.3",
      8.00, 0.15, 0.30, 2.0, 0.99, 0.30, 0.98, 8.0, 35
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

  TNamed promptKinematicsInfo(
    "prompt_kinematics",
    usePrimaryVertexKinematicsForPrompt
      ? "daughter momenta at independent transverse PCAs to the configured beam axis"
      : "secondary pair-PCA daughter momenta");
  promptKinematicsInfo.Write();

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

  // Prompt-meson selections are intentionally separate from the V0 cuts.
  // These candidates are treated as primary within the present DCA resolution.
  const std::vector<PromptSelection> promptSelections = {
    {
      "primary_dz_2p0",
      "primary momenta, pT>0.20, nTPC>=20, quality<20, |DCAxy|<3, "
      "|primary PCA z|<20, |Delta primary PCA z|<2.0; DCAz/PID/prompt-pair-DCA disabled",
      0.20, 20, 20.0, 3.0, -1.0, 20.0, 2.0, -1.0
    },
    {
      "primary_dz_1p0",
      "primary momenta, pT>0.20, nTPC>=20, quality<20, |DCAxy|<3, "
      "|primary PCA z|<20, |Delta primary PCA z|<1.0; DCAz/PID/prompt-pair-DCA disabled",
      0.20, 20, 20.0, 3.0, -1.0, 20.0, 1.0, -1.0
    },
    {
      "primary_dz_0p5",
      "primary momenta, pT>0.20, nTPC>=20, quality<20, |DCAxy|<3, "
      "|primary PCA z|<20, |Delta primary PCA z|<0.5; DCAz/PID/prompt-pair-DCA disabled",
      0.20, 20, 20.0, 3.0, -1.0, 20.0, 0.5, -1.0
    },
    {
      "primary_dz_0p2",
      "primary momenta, pT>0.20, nTPC>=20, quality<20, |DCAxy|<3, "
      "|primary PCA z|<20, |Delta primary PCA z|<0.2; DCAz/PID/prompt-pair-DCA disabled",
      0.20, 20, 20.0, 3.0, -1.0, 20.0, 0.2, -1.0
    },
    {
      "primary_dz_0p5_track",
      "primary momenta, pT>0.30, nTPC>=30, quality<15, |DCAxy|<1, "
      "|primary PCA z|<20, |Delta primary PCA z|<0.5",
      0.30, 30, 15.0, 1.0, -1.0, 20.0, 0.5, -1.0
    },
    {
      "primary_dz_0p5_dcaz2",
      "same as primary_dz_0p5_track plus |DCAz to fixed z=0|<2; diagnostic only",
      0.30, 30, 15.0, 1.0, 2.0, 20.0, 0.5, -1.0
    },
    {
      "primary_dz_0p5_pid",
      "same as primary_dz_0p5_track plus pion dE/dx<300 and kaon dE/dx>300",
      0.30, 30, 15.0, 1.0, -1.0, 20.0, 0.5, -1.0,
      -1.0, 300.0, 300.0, -1.0
    }
  };

  std::map<std::string, std::map<std::string, HistSet>> promptHistograms;

  TDirectory* promptTopDir = output->mkdir("promptMesons");
  for (const auto& selection : promptSelections)
  {
    TDirectory* selectionDir = promptTopDir->mkdir(selection.name.c_str());
    selectionDir->cd();

    TNamed selectionInfo("selection", selection.description.c_str());
    selectionInfo.Write();

    Selection bookingSelection{
      selection.name, selection.description,
      0., 0., 0., 0., 0., 0., 0., 0., 0
    };

    for (const auto category : {
           ChargeCategory::Unlike,
           ChargeCategory::Like,
           ChargeCategory::PlusPlus,
           ChargeCategory::MinusMinus})
    {
      const std::string categoryName = chargeName(category);
      TDirectory* chargeDir = selectionDir->mkdir(categoryName.c_str());
      promptHistograms[selection.name][categoryName] =
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
    << ", prompt kinematics = "
    << (usePrimaryVertexKinematicsForPrompt
          ? "primary vertex"
          : "secondary pair PCA")
    << std::endl;

  for (Long64_t entry = 0;
       entry < entriesToProcess;
       ++entry)
  {
    chain.GetEntry(entry);

    const double qualityScale =
      ScaleQualityBy10 ? 10.0 : 1.0;

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
    const bool passV0DeltaPhi =
      deltaPhi >= 0.8 - 0.4 * (v0_pt < 2.0 ? v0_pt : 2.0);

    const auto categories =
      chargeCategories(charge1, charge2);

    if (categories.empty())
    {
      continue;
    }

    // Secondary pair-PCA kinematics: always used for K0S/Lambda.
    const double pt1 = std::hypot(px1, py1);
    const double pt2 = std::hypot(px2, py2);

    const double daughterPtMin =
      std::min(pt1, pt2);

    // Prompt kinematics: selectable between the primary vertex and pair PCA.
    const double promptPx1 =
      usePrimaryVertexKinematicsForPrompt ? primary_px1 : px1;
    const double promptPy1 =
      usePrimaryVertexKinematicsForPrompt ? primary_py1 : py1;
    const double promptPz1 =
      usePrimaryVertexKinematicsForPrompt ? primary_pz1 : pz1;

    const double promptPx2 =
      usePrimaryVertexKinematicsForPrompt ? primary_px2 : px2;
    const double promptPy2 =
      usePrimaryVertexKinematicsForPrompt ? primary_py2 : py2;
    const double promptPz2 =
      usePrimaryVertexKinematicsForPrompt ? primary_pz2 : pz2;

    const bool promptKinematicsFinite =
      primary_pca_valid != 0 &&
      std::isfinite(primary_pca1_z) &&
      std::isfinite(primary_pca2_z) &&
      std::isfinite(promptPx1) &&
      std::isfinite(promptPy1) &&
      std::isfinite(promptPz1) &&
      std::isfinite(promptPx2) &&
      std::isfinite(promptPy2) &&
      std::isfinite(promptPz2);

    const double promptPt1 =
      std::hypot(promptPx1, promptPy1);
    const double promptPt2 =
      std::hypot(promptPx2, promptPy2);
    const double promptDaughterPtMin =
      std::min(promptPt1, promptPt2);

    const double promptPairPt =
      std::hypot(promptPx1 + promptPx2,
                 promptPy1 + promptPy2);

    const double promptMassPhi =
      invariantMass(
        promptPx1, promptPy1, promptPz1, kKaonMass,
        promptPx2, promptPy2, promptPz2, kKaonMass);

    double promptMassD0 = -1.0;
    double promptMassAntiD0 = -1.0;

    if (charge1 > 0)
    {
      // D0: K- pi+; anti-D0: K+ pi-.
      promptMassD0 =
        invariantMass(
          promptPx2, promptPy2, promptPz2, kKaonMass,
          promptPx1, promptPy1, promptPz1, kPionMass);

      promptMassAntiD0 =
        invariantMass(
          promptPx1, promptPy1, promptPz1, kKaonMass,
          promptPx2, promptPy2, promptPz2, kPionMass);
    }
    else
    {
      promptMassD0 =
        invariantMass(
          promptPx1, promptPy1, promptPz1, kKaonMass,
          promptPx2, promptPy2, promptPz2, kPionMass);

      promptMassAntiD0 =
        invariantMass(
          promptPx2, promptPy2, promptPz2, kKaonMass,
          promptPx1, promptPy1, promptPz1, kPionMass);
    }

    const double absDeltaPcaZ =
      std::abs(pca1_z - pca2_z);

    // Secondary V0 flight is measured from the configured beam position,
    // not from the detector origin.
    const double flightX = pca_x - beamX;
    const double flightY = pca_y - beamY;
    const double flightZ = pca_z - beamZ;

    const double decayRadius =
      std::hypot(flightX, flightY);

    const double absAlpha =
      std::abs(alpha);

    const double absPairDCA =
      std::abs(pairDCA);

    const double maxAbsTrackDcaXY =
      std::max(std::abs(dca_xy1), std::abs(dca_xy2));
    const double maxAbsTrackDcaZ =
      std::max(std::abs(dca_z1), std::abs(dca_z2));

    const double qualityMax =
      qualityScale * std::max(quality1, quality2);

    const int npointsMin =
      std::min<int>(npoints1, npoints2);

    const unsigned int ntpcClustersMin =
      std::min(ntpc_clusters1, ntpc_clusters2);

    const double maxAbsPrimaryPcaZ =
      std::max(
        std::abs(primary_pca1_z - beamZ),
        std::abs(primary_pca2_z - beamZ));

    const double absDeltaPrimaryPcaZ =
      std::isfinite(primary_pca_dz)
        ? std::abs(primary_pca_dz)
        : std::abs(primary_pca1_z - primary_pca2_z);

    const double absPromptPairDCA =
      (prompt_pca_valid != 0 && std::isfinite(prompt_pairDCA))
        ? std::abs(prompt_pairDCA)
        : std::numeric_limits<double>::quiet_NaN();

    const double v0Momentum =
      std::sqrt(
        v0_px * v0_px +
        v0_py * v0_py +
        v0_pz * v0_pz);

    const double flightLength =
      std::sqrt(
        flightX * flightX +
        flightY * flightY +
        flightZ * flightZ);

    const double dira =
      (v0Momentum > 0. && flightLength > 0.)
        ? (v0_px * flightX +
           v0_py * flightY +
           v0_pz * flightZ) /
          (v0Momentum * flightLength)
        : -2.;

    const bool unlikeSign =
      charge1 * charge2 < 0;

    const bool likeSign =
      charge1 * charge2 > 0;

    double positivePt = -1.;
    double negativePt = -1.;
    double positiveDedx = -1.;
    double negativeDedx = -1.;

    if (unlikeSign)
    {
      if (charge1 > 0)
      {
        positivePt = pt1;
        negativePt = pt2;
        positiveDedx = dedx_1;
        negativeDedx = dedx_2;
      }
      else
      {
        positivePt = pt2;
        negativePt = pt1;
        positiveDedx = dedx_2;
        negativeDedx = dedx_1;
      }
    }

    const int effectiveNpoints1 = npoints1;

    const int effectiveNpoints2 = npoints2;

    const bool exactCut1 =
      passV0DeltaPhi &&
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
      passV0DeltaPhi &&
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
      passV0DeltaPhi &&
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

      const bool kshortDedxPass =
        passesDedx(dedx_1, selection.pionDedxMin, selection.pionDedxMax) &&
        passesDedx(dedx_2, selection.pionDedxMin, selection.pionDedxMax);

      const bool lambdaDedxPass =
        unlikeSign &&
        passesDedx(positiveDedx,
                   selection.protonDedxMin,
                   selection.protonDedxMax) &&
        passesDedx(negativeDedx,
                   selection.pionDedxMin,
                   selection.pionDedxMax);

      const bool antiLambdaDedxPass =
        unlikeSign &&
        passesDedx(positiveDedx,
                   selection.pionDedxMin,
                   selection.pionDedxMax) &&
        passesDedx(negativeDedx,
                   selection.protonDedxMin,
                   selection.protonDedxMax);

      //if (!passV0DeltaPhi)
      //{
      //  continue;
      //}

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

        if (kshortDedxPass && passV0DeltaPhi)
        {
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
        }

        h.h_armenteros_podolanski->Fill(
          alpha, qT);

        h.h_pair_dca_vs_delta_pca_z->Fill(
          absDeltaPcaZ, absPairDCA);

        // Lambda and anti-Lambda are physically meaningful only
        // for unlike-sign pairs.
        if (category == ChargeCategory::Unlike)
        {
          if (lambdaProtonPtPass && lambdaDedxPass && passV0DeltaPhi)
          {
            h.h_lambda_mass->Fill(mass_Lambda);
            h.h_lambda_mass_vs_v0pt->Fill(
              v0_pt, mass_Lambda);

            h.h3_lambda->Fill(
              v0_pt, mass_Lambda, absPairDCA);
          }

          if (antiLambdaProtonPtPass && antiLambdaDedxPass && passV0DeltaPhi)
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

    // Prompt phi and D0/anti-D0 QA.
    // Physical prompt candidates are unlike-sign. The present reconstruction
    // retains same-sign rows only through the K0S background path, so those
    // rows are not an unbiased prompt like-sign sample and are not filled here.
    for (const auto& promptSelection : promptSelections)
    {
      const bool phiDedxPass =
        passesDedx(dedx_1,
                   promptSelection.kaonDedxMin,
                   promptSelection.kaonDedxMax) &&
        passesDedx(dedx_2,
                   promptSelection.kaonDedxMin,
                   promptSelection.kaonDedxMax);

      const bool d0DedxPass =
        unlikeSign &&
        passesDedx(positiveDedx,
                   promptSelection.pionDedxMin,
                   promptSelection.pionDedxMax) &&
        passesDedx(negativeDedx,
                   promptSelection.kaonDedxMin,
                   promptSelection.kaonDedxMax);

      const bool antiD0DedxPass =
        unlikeSign &&
        passesDedx(positiveDedx,
                   promptSelection.kaonDedxMin,
                   promptSelection.kaonDedxMax) &&
        passesDedx(negativeDedx,
                   promptSelection.pionDedxMin,
                   promptSelection.pionDedxMax);

      if (!promptKinematicsFinite ||
          !passesPromptSelection(
            promptSelection,
            promptDaughterPtMin,
            ntpcClustersMin,
            qualityMax,
            maxAbsTrackDcaXY,
            maxAbsTrackDcaZ,
            maxAbsPrimaryPcaZ,
            absDeltaPrimaryPcaZ,
            absPromptPairDCA))
      {
        continue;
      }

      for (const auto category : categories)
      {
        if (category != ChargeCategory::Unlike)
        {
          continue;
        }

        HistSet& h =
          promptHistograms
            .at(promptSelection.name)
            .at(chargeName(category));

        if (phiDedxPass)
        {
          h.h_phi_mass->Fill(promptMassPhi);
          h.h_phi_mass_vs_v0pt->Fill(
            promptPairPt, promptMassPhi);

          if (std::isfinite(absPromptPairDCA))
          {
            h.h3_phi->Fill(
              promptPairPt, promptMassPhi, absPromptPairDCA);
          }
        }

        h.h_primary_pca_dz->Fill(absDeltaPrimaryPcaZ);
        h.h_primary_pca_z1_vs_z2->Fill(
          primary_pca1_z, primary_pca2_z);
        h.h_max_dca_xy_vs_primary_pca_dz->Fill(
          absDeltaPrimaryPcaZ, maxAbsTrackDcaXY);

        if (std::isfinite(absPromptPairDCA))
        {
          h.h_pair_dca_vs_delta_pca_z->Fill(
            absDeltaPrimaryPcaZ, absPromptPairDCA);
        }

        if (category == ChargeCategory::Unlike)
        {
          if (d0DedxPass)
          {
            h.h_d0_mass->Fill(promptMassD0);
            h.h_d0_mass_vs_v0pt->Fill(
              promptPairPt, promptMassD0);
            if (std::isfinite(absPromptPairDCA))
            {
              h.h3_d0->Fill(
                promptPairPt, promptMassD0, absPromptPairDCA);
            }
          }

          if (antiD0DedxPass)
          {
            h.h_antid0_mass->Fill(promptMassAntiD0);
            h.h_antid0_mass_vs_v0pt->Fill(
              promptPairPt, promptMassAntiD0);
            if (std::isfinite(absPromptPairDCA))
            {
              h.h3_antid0->Fill(
                promptPairPt, promptMassAntiD0, absPromptPairDCA);
            }
          }
        }
      }
    }

    const bool exactKshortDedxPass =
      passesDedx(dedx_1, -1.0, 400.0) &&
      passesDedx(dedx_2, -1.0, 400.0);

    const bool exactLambdaDedxPass =
      unlikeSign && passesDedx(negativeDedx, -1.0, 400.0);

    const bool exactAntiLambdaDedxPass =
      unlikeSign && passesDedx(positiveDedx, -1.0, 400.0);

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

        if (exactKshortDedxPass)
        {
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
        }

        h.h_armenteros_podolanski->Fill(alpha, qT);
        h.h_pair_dca_vs_delta_pca_z->Fill(
          absDeltaPcaZ, absPairDCA);

        if (category == ChargeCategory::Unlike)
        {
          if (lambdaProtonPtPass && exactLambdaDedxPass)
          {
            h.h_lambda_mass->Fill(mass_Lambda);
            h.h_lambda_mass_vs_v0pt->Fill(v0_pt, mass_Lambda);
            h.h3_lambda->Fill(v0_pt, mass_Lambda, absPairDCA);
          }

          if (antiLambdaProtonPtPass && exactAntiLambdaDedxPass)
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
