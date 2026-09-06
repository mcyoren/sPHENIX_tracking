// MakeSigmaCascadeHistograms_extended.C
//
// One-pass QA/cut scan for TpcHadronReco:
//   sigma1385Tree : Sigma(1385) -> Lambda pi
//   cascadeTree   : Xi -> Lambda pi, Omega -> Lambda K
//
// This version KEEPS the original cut scan and original histograms, and adds:
//
//   * bachelor dE/dx QA;
//   * ALEPH-style bachelor PID categories:
//         pid_all
//         pid_loose   |N_sigma| < 3
//         pid_normal  |N_sigma| < 2
//         pid_tight   |N_sigma| < 1
//     Sigma/Xi use the pion hypothesis; Omega uses the kaon hypothesis.
//
//   * Lambda-bachelor angular QA:
//         opening angle in 3D
//         signed Delta-phi
//         |Delta-eta|
//         |Delta-phi|/(|Delta-eta|+epsilon)
//
//     The signed Delta-phi is charge oriented:
//       positive bachelor:  phi_positive - phi_Lambda
//       negative bachelor:  phi_Lambda   - phi_negative
//
//     This is the direct cascade analogue of the signed opening requirement
//     used in the V0 QA.
//
//   * within EVERY original cut directory and EVERY PID category:
//         mass
//         mass after loose/normal/tight signed-angle cuts
//         mass after phi-dominated opening
//         mass vs signed Delta-phi
//         mass vs opening angle
//         mass vs |Delta-eta|
//         mass vs |Delta-phi|/(|Delta-eta|+eps)
//         mass vs bachelor p
//         N_sigma(target) vs bachelor p
//         ln(dE/dx) vs bachelor p
//         TH3F: mass vs mother pT vs signed Delta-phi
//
// The original producer selections are already present in the trees, so this
// macro only tightens / partitions what was written.
//
// Example:
// root -l -b -q \
// 'MakeSigmaCascadeHistograms_extended.C("/path","v0_pp_*.root","qa","sigma_cascade_extended.root")'
//

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>

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
constexpr double kLambdaMass = 1.115683;
constexpr double kPionMass   = 0.13957039;
constexpr double kKaonMass   = 0.493677;
constexpr double kTwoPi      = 2.0 * TMath::Pi();

constexpr UInt_t kXi    = 1U << 0;
constexpr UInt_t kOmega = 1U << 1;

// ALEPH/ALICE central dE/dx curve used in the previous PID study.
constexpr double kAlephP1 = 0.1819;
constexpr double kAlephP2 = 10.8175;
constexpr double kAlephP3 = 1.1674e-5;
constexpr double kAlephP4 = 3.0199;
constexpr double kAlephP5 = 5.0;

// Universal ln(dE/dx) width.
constexpr double kWidthA = 0.07690;
constexpr double kWidthB = 0.16914;
constexpr double kWidthC = 0.40000;

bool has(TChain& t, const char* b)
{
  return t.GetBranch(b) != nullptr;
}

double wrapPhi(double x)
{
  while (x >= TMath::Pi()) x -= kTwoPi;
  while (x < -TMath::Pi()) x += kTwoPi;
  return x;
}

double mag3(double px, double py, double pz)
{
  return std::sqrt(px*px + py*py + pz*pz);
}

double eta(double px, double py, double pz)
{
  const double pt = std::hypot(px, py);
  if (!(pt > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  return std::asinh(pz/pt);
}

double openingAngle(double px1, double py1, double pz1,
                    double px2, double py2, double pz2)
{
  const double p1 = mag3(px1, py1, pz1);
  const double p2 = mag3(px2, py2, pz2);

  if (!(p1 > 0.0) || !(p2 > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  double cosine =
    (px1*px2 + py1*py2 + pz1*pz2)/(p1*p2);

  cosine = std::clamp(cosine, -1.0, 1.0);
  return std::acos(cosine);
}

double alephExpectedDedx(double betaGamma)
{
  if (!(betaGamma > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  const double beta =
    betaGamma/std::sqrt(1.0 + betaGamma*betaGamma);

  const double betaP4 =
    std::pow(beta, kAlephP4);

  const double argument =
    kAlephP3 + std::pow(betaGamma, -kAlephP5);

  if (!(betaP4 > 0.0) || !(argument > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  const double value =
    (kAlephP1/betaP4) *
    (kAlephP2 - betaP4 - std::log(argument));

  return value > 0.0
    ? value
    : std::numeric_limits<double>::quiet_NaN();
}

double universalWidth(double expectedDedx)
{
  if (!(expectedDedx > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  return std::sqrt(
    kWidthA*kWidthA +
    kWidthB*kWidthB/std::pow(expectedDedx, kWidthC));
}

// The pion kernel-smoothed width points were not tabulated in the presentation,
// so for pion we keep the universal width.  For kaons use the published
// momentum-dependent width-ratio correction, frozen outside its fit range.
double kaonWidthRatio(double momentum)
{
  constexpr double pLow  = 0.23;
  constexpr double pHigh = 0.53;

  const double p =
    std::clamp(momentum, pLow, pHigh);

  return
    0.910 +
    0.628*(std::log(p) - std::log(pLow));
}

double targetNSigma(double momentum,
                    double dedx,
                    double mass,
                    bool useKaonWidthCorrection,
                    double dedxScale)
{
  if (!(momentum > 0.0) ||
      !(dedx > 0.0) ||
      !(mass > 0.0) ||
      !(dedxScale > 0.0))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double expected =
    alephExpectedDedx(momentum/mass);

  if (!(expected > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  double sigma =
    universalWidth(expected);

  if (useKaonWidthCorrection)
    sigma *= kaonWidthRatio(momentum);

  if (!(sigma > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  return
    (std::log(dedx*dedxScale) - std::log(expected))/sigma;
}

// Charge-oriented signed opening:
//
// bachelor q > 0 : phi+ - phi_Lambda
// bachelor q < 0 : phi_Lambda - phi-
//
// This keeps the expected bending/opening direction aligned between particle
// and antiparticle candidates.
double signedLambdaBachelorDeltaPhi(double lambdaPx,
                                    double lambdaPy,
                                    double bachelorPx,
                                    double bachelorPy,
                                    int bachelorCharge)
{
  const double phiLambda =
    std::atan2(lambdaPy, lambdaPx);

  const double phiBachelor =
    std::atan2(bachelorPy, bachelorPx);

  return bachelorCharge > 0
    ? wrapPhi(phiBachelor - phiLambda)
    : wrapPhi(phiLambda - phiBachelor);
}

// A deliberately mild, tunable analogue of the K0S pT-dependent signed-angle
// cut.  "atZero" is the threshold at mother pT=0 and it decreases linearly to
// zero at mother pT=zeroAtPt.
//
// These angle subsets are QA only: the original cut folders remain untouched.
double signedAngleThreshold(double motherPt,
                            double atZero,
                            double zeroAtPt)
{
  if (!(zeroAtPt > 0.0))
    return std::max(0.0, atZero);

  const double fraction =
    std::clamp(motherPt/zeroAtPt, 0.0, 1.0);

  return
    std::max(0.0, atZero*(1.0 - fraction));
}

struct SigmaCut
{
  std::string name, desc;
  double dmLambda, maxLambdaDca, minLambdaDira, minLambdaR;
  double minBachelorPt;
  int minBachelorN;
  double maxBachelorDcaXY;
};

struct CascadeCut
{
  std::string name, desc;
  double dmLambda, maxLambdaDca, minBachelorPt;
  int minBachelorN;
  double maxCascadeDca, minCascadeR, minCascadeDira;
  double minLambdaFlight, minLambdaDira, minBachelorDcaXY;
};

bool passSigma(const SigmaCut& c,
               double ml,
               double ldca,
               double ldira,
               double lr,
               double bpt,
               unsigned int bn,
               double bdca)
{
  if (c.dmLambda >= 0 &&
      (!std::isfinite(ml) ||
       std::abs(ml-kLambdaMass) > c.dmLambda))
    return false;

  if (c.maxLambdaDca >= 0 &&
      (!std::isfinite(ldca) ||
       std::abs(ldca) > c.maxLambdaDca))
    return false;

  if (c.minLambdaDira >= -1 &&
      (!std::isfinite(ldira) ||
       ldira < c.minLambdaDira))
    return false;

  if (c.minLambdaR >= 0 &&
      (!std::isfinite(lr) ||
       lr < c.minLambdaR))
    return false;

  if (c.minBachelorPt >= 0 &&
      (!std::isfinite(bpt) ||
       bpt < c.minBachelorPt))
    return false;

  if (c.minBachelorN > 0 &&
      static_cast<int>(bn) < c.minBachelorN)
    return false;

  if (c.maxBachelorDcaXY >= 0 &&
      (!std::isfinite(bdca) ||
       std::abs(bdca) > c.maxBachelorDcaXY))
    return false;

  return true;
}

bool passCascade(const CascadeCut& c,
                 double ml,
                 double ldca,
                 double bpt,
                 unsigned int bn,
                 double cdca,
                 double cr,
                 double cdira,
                 double lflight,
                 double ldira,
                 double bdca)
{
  if (c.dmLambda >= 0 &&
      (!std::isfinite(ml) ||
       std::abs(ml-kLambdaMass) > c.dmLambda))
    return false;

  if (c.maxLambdaDca >= 0 &&
      (!std::isfinite(ldca) ||
       std::abs(ldca) > c.maxLambdaDca))
    return false;

  if (c.minBachelorPt >= 0 &&
      (!std::isfinite(bpt) ||
       bpt < c.minBachelorPt))
    return false;

  if (c.minBachelorN > 0 &&
      static_cast<int>(bn) < c.minBachelorN)
    return false;

  if (c.maxCascadeDca >= 0 &&
      (!std::isfinite(cdca) ||
       std::abs(cdca) > c.maxCascadeDca))
    return false;

  if (c.minCascadeR >= 0 &&
      (!std::isfinite(cr) ||
       cr < c.minCascadeR))
    return false;

  if (c.minCascadeDira >= -1 &&
      (!std::isfinite(cdira) ||
       cdira < c.minCascadeDira))
    return false;

  if (c.minLambdaFlight >= 0 &&
      (!std::isfinite(lflight) ||
       lflight < c.minLambdaFlight))
    return false;

  if (c.minLambdaDira >= -1 &&
      (!std::isfinite(ldira) ||
       ldira < c.minLambdaDira))
    return false;

  if (c.minBachelorDcaXY >= 0 &&
      (!std::isfinite(bdca) ||
       std::abs(bdca) < c.minBachelorDcaXY))
    return false;

  return true;
}

// ============================================================================
// Original histogram blocks
// ============================================================================

struct SigmaH
{
  TH1F *mass{}, *lamP{}, *lamM{}, *alamP{}, *alamM{};

  TH2F
    *pt{},
    *lambdaMass{},
    *lambdaR{},
    *lambdaDca{},
    *lambdaDira{},
    *bachelorPt{},
    *bachelorDca{};
};

SigmaH bookSigma(TDirectory* d,
                 const std::string& tag)
{
  d->cd();

  SigmaH h;

  h.mass = new TH1F(
    "h_mass_sigma1385",
    ("#Sigma(1385) #rightarrow #Lambda#pi ["+tag+
     "];m_{#Lambda#pi} [GeV/c^{2}];candidates").c_str(),
    300, 1.25, 1.55);

  h.lamP = new TH1F(
    "h_mass_LambdaPiPlus",
    "#Lambda#pi^{+};m [GeV/c^{2}];candidates",
    300, 1.25, 1.55);

  h.lamM = new TH1F(
    "h_mass_LambdaPiMinus",
    "#Lambda#pi^{-};m [GeV/c^{2}];candidates",
    300, 1.25, 1.55);

  h.alamP = new TH1F(
    "h_mass_AntiLambdaPiPlus",
    "#bar{#Lambda}#pi^{+};m [GeV/c^{2}];candidates",
    300, 1.25, 1.55);

  h.alamM = new TH1F(
    "h_mass_AntiLambdaPiMinus",
    "#bar{#Lambda}#pi^{-};m [GeV/c^{2}];candidates",
    300, 1.25, 1.55);

  h.pt = new TH2F(
    "h_mass_vs_pt",
    "mass vs p_{T};p_{T}^{#Lambda#pi} [GeV/c];m [GeV/c^{2}]",
    60, 0, 6,
    300, 1.25, 1.55);

  h.lambdaMass = new TH2F(
    "h_mass_vs_lambda_mass",
    "mass vs #Lambda mass;m_{#Lambda};m_{#Lambda#pi}",
    160, 1.075, 1.155,
    300, 1.25, 1.55);

  h.lambdaR = new TH2F(
    "h_mass_vs_lambda_decayR",
    "mass vs #Lambda decay radius;R_{#Lambda} [cm];m",
    100, 0, 50,
    300, 1.25, 1.55);

  h.lambdaDca = new TH2F(
    "h_mass_vs_lambda_pairDCA",
    "mass vs #Lambda pair DCA;DCA_{p#pi} [cm];m",
    100, 0, 3,
    300, 1.25, 1.55);

  h.lambdaDira = new TH2F(
    "h_mass_vs_lambda_DIRA",
    "mass vs #Lambda DIRA;DIRA_{#Lambda};m",
    120, -0.2, 1.0,
    300, 1.25, 1.55);

  h.bachelorPt = new TH2F(
    "h_mass_vs_bachelor_pt",
    "mass vs bachelor p_{T};p_{T}^{#pi};m",
    60, 0, 3,
    300, 1.25, 1.55);

  h.bachelorDca = new TH2F(
    "h_mass_vs_bachelor_dcaXY",
    "mass vs bachelor |DCA_{xy}|;|DCA_{xy}| [cm];m",
    100, 0, 3,
    300, 1.25, 1.55);

  return h;
}

struct CascadeH
{
  TH1F *xi{}, *xim{}, *xip{}, *om{}, *omm{}, *omp{};

  TH2F
    *xiPt{},
    *omPt{},
    *xiDca{},
    *omDca{},
    *xiR{},
    *omR{},
    *xiDira{},
    *omDira{},
    *xiLFlight{},
    *omLFlight{},
    *xiLDira{},
    *omLDira{},
    *xiBDca{},
    *omBDca{},
    *xiLmass{},
    *omLmass{};
};

CascadeH bookCascade(TDirectory* d)
{
  d->cd();

  CascadeH h;

  h.xi = new TH1F(
    "h_mass_xi",
    "#Xi^{#mp};m_{#Lambda#pi} [GeV/c^{2}];candidates",
    300, 1.25, 1.40);

  h.xim = new TH1F(
    "h_mass_xi_minus",
    "#Xi^{-};m_{#Lambda#pi^{-}} [GeV/c^{2}];candidates",
    300, 1.25, 1.40);

  h.xip = new TH1F(
    "h_mass_xi_plus",
    "#bar{#Xi}^{+};m_{#bar{#Lambda}#pi^{+}} [GeV/c^{2}];candidates",
    300, 1.25, 1.40);

  h.om = new TH1F(
    "h_mass_omega",
    "#Omega^{#mp};m_{#Lambda K} [GeV/c^{2}];candidates",
    300, 1.60, 1.75);

  h.omm = new TH1F(
    "h_mass_omega_minus",
    "#Omega^{-};m_{#Lambda K^{-}} [GeV/c^{2}];candidates",
    300, 1.60, 1.75);

  h.omp = new TH1F(
    "h_mass_omega_plus",
    "#bar{#Omega}^{+};m_{#bar{#Lambda}K^{+}} [GeV/c^{2}];candidates",
    300, 1.60, 1.75);

  h.xiPt = new TH2F(
    "h_xi_mass_vs_pt",
    ";p_{T}^{#Xi};m_{#Xi}",
    60, 0, 6,
    300, 1.25, 1.40);

  h.omPt = new TH2F(
    "h_omega_mass_vs_pt",
    ";p_{T}^{#Omega};m_{#Omega}",
    60, 0, 6,
    300, 1.60, 1.75);

  h.xiDca = new TH2F(
    "h_xi_mass_vs_cascadeDCA",
    ";DCA_{#Lambda,bach} [cm];m_{#Xi}",
    100, 0, 2,
    300, 1.25, 1.40);

  h.omDca = new TH2F(
    "h_omega_mass_vs_cascadeDCA",
    ";DCA_{#Lambda,bach} [cm];m_{#Omega}",
    100, 0, 2,
    300, 1.60, 1.75);

  h.xiR = new TH2F(
    "h_xi_mass_vs_decayR",
    ";R_{cascade} [cm];m_{#Xi}",
    100, 0, 50,
    300, 1.25, 1.40);

  h.omR = new TH2F(
    "h_omega_mass_vs_decayR",
    ";R_{cascade} [cm];m_{#Omega}",
    100, 0, 50,
    300, 1.60, 1.75);

  h.xiDira = new TH2F(
    "h_xi_mass_vs_DIRA",
    ";DIRA_{cascade};m_{#Xi}",
    100, 0.70, 1,
    300, 1.25, 1.40);

  h.omDira = new TH2F(
    "h_omega_mass_vs_DIRA",
    ";DIRA_{cascade};m_{#Omega}",
    100, 0.70, 1,
    300, 1.60, 1.75);

  h.xiLFlight = new TH2F(
    "h_xi_mass_vs_lambdaFlight",
    ";L_{#Lambda}^{cascade} [cm];m_{#Xi}",
    100, 0, 50,
    300, 1.25, 1.40);

  h.omLFlight = new TH2F(
    "h_omega_mass_vs_lambdaFlight",
    ";L_{#Lambda}^{cascade} [cm];m_{#Omega}",
    100, 0, 50,
    300, 1.60, 1.75);

  h.xiLDira = new TH2F(
    "h_xi_mass_vs_lambdaDIRA",
    ";DIRA_{#Lambda}^{cascade};m_{#Xi}",
    100, 0.70, 1,
    300, 1.25, 1.40);

  h.omLDira = new TH2F(
    "h_omega_mass_vs_lambdaDIRA",
    ";DIRA_{#Lambda}^{cascade};m_{#Omega}",
    100, 0.70, 1,
    300, 1.60, 1.75);

  h.xiBDca = new TH2F(
    "h_xi_mass_vs_bachelorDCAxy",
    ";|DCA_{xy}^{bach}| [cm];m_{#Xi}",
    100, 0, 5,
    300, 1.25, 1.40);

  h.omBDca = new TH2F(
    "h_omega_mass_vs_bachelorDCAxy",
    ";|DCA_{xy}^{bach}| [cm];m_{#Omega}",
    100, 0, 5,
    300, 1.60, 1.75);

  h.xiLmass = new TH2F(
    "h_xi_mass_vs_lambdaMass",
    ";m_{#Lambda};m_{#Xi}",
    160, 1.075, 1.155,
    300, 1.25, 1.40);

  h.omLmass = new TH2F(
    "h_omega_mass_vs_lambdaMass",
    ";m_{#Lambda};m_{#Omega}",
    160, 1.075, 1.155,
    300, 1.60, 1.75);

  return h;
}

// ============================================================================
// New PID + angular QA block
// ============================================================================

struct AngularPidH
{
  TH1F
    *mass{},
    *massAngleLoose{},
    *massAngleNormal{},
    *massAngleTight{},
    *massPhiDominated{},
    *massAngleNormalPhiDominated{};

  TH2F
    *massVsSignedDphi{},
    *massVsOpening{},
    *massVsAbsDeta{},
    *massVsPhiEtaRatio{},
    *massVsBachelorP{},
    *nSigmaVsBachelorP{},
    *lnDedxVsBachelorP{};

  TH3F
    *massVsPtVsSignedDphi{};
};

AngularPidH bookAngularPid(TDirectory* d,
                           const std::string& title,
                           double massMin,
                           double massMax,
                           const std::string& bachelorName,
                           const std::string& pidName)
{
  d->cd();

  AngularPidH h;

  h.mass = new TH1F(
    "h_mass",
    (title+" ["+pidName+"];mass [GeV/c^{2}];candidates").c_str(),
    300, massMin, massMax);

  h.massAngleLoose = new TH1F(
    "h_mass_angle_loose",
    (title+" + loose signed-angle cut;mass [GeV/c^{2}];candidates").c_str(),
    300, massMin, massMax);

  h.massAngleNormal = new TH1F(
    "h_mass_angle_normal",
    (title+" + normal signed-angle cut;mass [GeV/c^{2}];candidates").c_str(),
    300, massMin, massMax);

  h.massAngleTight = new TH1F(
    "h_mass_angle_tight",
    (title+" + tight signed-angle cut;mass [GeV/c^{2}];candidates").c_str(),
    300, massMin, massMax);

  h.massPhiDominated = new TH1F(
    "h_mass_phi_dominated",
    (title+" + |#Delta#phi|>|#Delta#eta|;mass [GeV/c^{2}];candidates").c_str(),
    300, massMin, massMax);

  h.massAngleNormalPhiDominated = new TH1F(
    "h_mass_angle_normal_phi_dominated",
    (title+" + normal signed-angle + #phi-dominated opening;mass [GeV/c^{2}];candidates").c_str(),
    300, massMin, massMax);

  h.massVsSignedDphi = new TH2F(
    "h_mass_vs_signed_lambda_bachelor_dphi",
    (title+
     ";signed #Delta#phi_{#Lambda,bach} [rad];mass [GeV/c^{2}]").c_str(),
    144, -TMath::Pi(), TMath::Pi(),
    300, massMin, massMax);

  h.massVsOpening = new TH2F(
    "h_mass_vs_lambda_bachelor_opening",
    (title+
     ";3D opening angle(#Lambda,"+bachelorName+") [rad];mass [GeV/c^{2}]").c_str(),
    120, 0.0, TMath::Pi(),
    300, massMin, massMax);

  h.massVsAbsDeta = new TH2F(
    "h_mass_vs_abs_lambda_bachelor_deta",
    (title+
     ";|#Delta#eta_{#Lambda,bach}|;mass [GeV/c^{2}]").c_str(),
    120, 0.0, 2.0,
    300, massMin, massMax);

  h.massVsPhiEtaRatio = new TH2F(
    "h_mass_vs_abs_dphi_over_abs_deta",
    (title+
     ";|#Delta#phi|/(|#Delta#eta|+10^{-3});mass [GeV/c^{2}]").c_str(),
    120, 0.0, 6.0,
    300, massMin, massMax);

  h.massVsBachelorP = new TH2F(
    "h_mass_vs_bachelor_p",
    (title+
     ";p_{"+bachelorName+"} [GeV/c];mass [GeV/c^{2}]").c_str(),
    100, 0.0, 5.0,
    300, massMin, massMax);

  h.nSigmaVsBachelorP = new TH2F(
    "h_nsigma_target_vs_bachelor_p",
    (title+
     ";p_{"+bachelorName+"} [GeV/c];N#sigma_{"+pidName+"}").c_str(),
    100, 0.0, 5.0,
    160, -8.0, 8.0);

  h.lnDedxVsBachelorP = new TH2F(
    "h_lnDedx_vs_bachelor_p",
    (title+
     ";p_{"+bachelorName+"} [GeV/c];ln(dE/dx)").c_str(),
    100, 0.0, 5.0,
    180, 0.0, 6.0);

  // Requested compact 3D:
  //
  //   mass vs mother pT vs
  //   (phi_positive - phi_Lambda) for q_bach>0
  //   (phi_Lambda - phi_negative) for q_bach<0
  //
  // Moderate binning keeps the output manageable even though this is booked
  // for every original cut and every PID category.
  h.massVsPtVsSignedDphi = new TH3F(
    "h3_mass_vs_pt_vs_signed_lambda_bachelor_dphi",
    (title+
     ";mass [GeV/c^{2}];mother p_{T} [GeV/c];"
     "signed #Delta#phi_{#Lambda,bach} [rad]").c_str(),
    120, massMin, massMax,
    30, 0.0, 6.0,
    72, -TMath::Pi(), TMath::Pi());

  return h;
}

void fillAngularPid(AngularPidH& h,
                    double mass,
                    double motherPt,
                    double bachelorP,
                    double bachelorDedx,
                    double nSigmaTarget,
                    double signedDphi,
                    double opening,
                    double absDeta,
                    double phiEtaRatio,
                    double angleLooseAtZero,
                    double angleNormalAtZero,
                    double angleTightAtZero,
                    double angleZeroAtPt)
{
  h.mass->Fill(mass);

  h.massVsSignedDphi->Fill(
    signedDphi,
    mass);

  if (std::isfinite(opening))
    h.massVsOpening->Fill(opening, mass);

  if (std::isfinite(absDeta))
    h.massVsAbsDeta->Fill(absDeta, mass);

  if (std::isfinite(phiEtaRatio))
    h.massVsPhiEtaRatio->Fill(phiEtaRatio, mass);

  if (std::isfinite(bachelorP))
    h.massVsBachelorP->Fill(bachelorP, mass);

  if (std::isfinite(nSigmaTarget))
    h.nSigmaVsBachelorP->Fill(bachelorP, nSigmaTarget);

  if (bachelorDedx > 0.0)
    h.lnDedxVsBachelorP->Fill(
      bachelorP,
      std::log(bachelorDedx));

  h.massVsPtVsSignedDphi->Fill(
    mass,
    motherPt,
    signedDphi);

  const double looseThreshold =
    signedAngleThreshold(
      motherPt,
      angleLooseAtZero,
      angleZeroAtPt);

  const double normalThreshold =
    signedAngleThreshold(
      motherPt,
      angleNormalAtZero,
      angleZeroAtPt);

  const double tightThreshold =
    signedAngleThreshold(
      motherPt,
      angleTightAtZero,
      angleZeroAtPt);

  if (signedDphi >= looseThreshold)
    h.massAngleLoose->Fill(mass);

  if (signedDphi >= normalThreshold)
    h.massAngleNormal->Fill(mass);

  if (signedDphi >= tightThreshold)
    h.massAngleTight->Fill(mass);

  const bool phiDominated =
    std::isfinite(absDeta) &&
    std::abs(signedDphi) > absDeta;

  if (phiDominated)
    h.massPhiDominated->Fill(mass);

  if (phiDominated &&
      signedDphi >= normalThreshold)
  {
    h.massAngleNormalPhiDominated->Fill(mass);
  }
}

struct PidFamily
{
  AngularPidH all;
  AngularPidH loose;
  AngularPidH normal;
  AngularPidH tight;
};

PidFamily bookPidFamily(TDirectory* cutDirectory,
                        const std::string& speciesPrefix,
                        const std::string& title,
                        double massMin,
                        double massMax,
                        const std::string& bachelorName,
                        const std::string& targetPid)
{
  PidFamily family;

  TDirectory* allDir =
    cutDirectory->mkdir(
      (speciesPrefix+"_pid_all").c_str());

  TDirectory* looseDir =
    cutDirectory->mkdir(
      (speciesPrefix+"_pid_loose").c_str());

  TDirectory* normalDir =
    cutDirectory->mkdir(
      (speciesPrefix+"_pid_normal").c_str());

  TDirectory* tightDir =
    cutDirectory->mkdir(
      (speciesPrefix+"_pid_tight").c_str());

  allDir->cd();
  TNamed(
    "pid_selection",
    "no bachelor PID cut").Write();

  looseDir->cd();
  TNamed(
    "pid_selection",
    ("|Nsigma_"+targetPid+"| < loose threshold").c_str()).Write();

  normalDir->cd();
  TNamed(
    "pid_selection",
    ("|Nsigma_"+targetPid+"| < normal threshold").c_str()).Write();

  tightDir->cd();
  TNamed(
    "pid_selection",
    ("|Nsigma_"+targetPid+"| < tight threshold").c_str()).Write();

  family.all =
    bookAngularPid(
      allDir,
      title,
      massMin,
      massMax,
      bachelorName,
      "all");

  family.loose =
    bookAngularPid(
      looseDir,
      title,
      massMin,
      massMax,
      bachelorName,
      targetPid);

  family.normal =
    bookAngularPid(
      normalDir,
      title,
      massMin,
      massMax,
      bachelorName,
      targetPid);

  family.tight =
    bookAngularPid(
      tightDir,
      title,
      massMin,
      massMax,
      bachelorName,
      targetPid);

  return family;
}

void fillPidFamily(PidFamily& family,
                   double mass,
                   double motherPt,
                   double bachelorP,
                   double bachelorDedxScaled,
                   double nSigmaTarget,
                   double signedDphi,
                   double opening,
                   double absDeta,
                   double phiEtaRatio,
                   double pidLooseNSigma,
                   double pidNormalNSigma,
                   double pidTightNSigma,
                   double angleLooseAtZero,
                   double angleNormalAtZero,
                   double angleTightAtZero,
                   double angleZeroAtPt)
{
  fillAngularPid(
    family.all,
    mass,
    motherPt,
    bachelorP,
    bachelorDedxScaled,
    nSigmaTarget,
    signedDphi,
    opening,
    absDeta,
    phiEtaRatio,
    angleLooseAtZero,
    angleNormalAtZero,
    angleTightAtZero,
    angleZeroAtPt);

  if (!std::isfinite(nSigmaTarget))
    return;

  if (std::abs(nSigmaTarget) < pidLooseNSigma)
  {
    fillAngularPid(
      family.loose,
      mass,
      motherPt,
      bachelorP,
      bachelorDedxScaled,
      nSigmaTarget,
      signedDphi,
      opening,
      absDeta,
      phiEtaRatio,
      angleLooseAtZero,
      angleNormalAtZero,
      angleTightAtZero,
      angleZeroAtPt);
  }

  if (std::abs(nSigmaTarget) < pidNormalNSigma)
  {
    fillAngularPid(
      family.normal,
      mass,
      motherPt,
      bachelorP,
      bachelorDedxScaled,
      nSigmaTarget,
      signedDphi,
      opening,
      absDeta,
      phiEtaRatio,
      angleLooseAtZero,
      angleNormalAtZero,
      angleTightAtZero,
      angleZeroAtPt);
  }

  if (std::abs(nSigmaTarget) < pidTightNSigma)
  {
    fillAngularPid(
      family.tight,
      mass,
      motherPt,
      bachelorP,
      bachelorDedxScaled,
      nSigmaTarget,
      signedDphi,
      opening,
      absDeta,
      phiEtaRatio,
      angleLooseAtZero,
      angleNormalAtZero,
      angleTightAtZero,
      angleZeroAtPt);
  }
}

} // namespace

// ============================================================================

void MakeSigmaCascadeHistograms(
    const char* inputDir=".",
    const char* filePattern="*.root",
    const char* outputDir="output",
    const char* outputName="sigma_cascade_histograms_extended.root",
    const char* sigmaTreeName="sigma1385Tree",
    const char* cascadeTreeName="cascadeTree",
    const Long64_t maxEntries=-1,
    const double beamX=0.158,
    const double beamY=0.285,

    // dE/dx scale and three nested PID working points.
    const double dedxScale=0.015,
    const double pidLooseNSigma=3.0,
    const double pidNormalNSigma=2.0,
    const double pidTightNSigma=1.0,

    // Exploratory signed-angle working points.
    //
    // threshold(pt) = threshold_at_zero * (1 - min(pt,zeroAtPt)/zeroAtPt)
    //
    // These DO NOT alter the existing cumulative cut scan; they create
    // additional mass spectra inside each folder.
    const double angleLooseAtZero=0.10,
    const double angleNormalAtZero=0.20,
    const double angleTightAtZero=0.30,
    const double angleZeroAtPt=2.0)
{
  TH1::AddDirectory(kTRUE);

  const TString pattern =
    TString::Format("%s/%s", inputDir, filePattern);

  gSystem->mkdir(outputDir, kTRUE);

  const TString outpath =
    TString::Format("%s/%s", outputDir, outputName);

  std::unique_ptr<TFile> out(
    TFile::Open(outpath, "RECREATE"));

  if (!out || out->IsZombie())
  {
    std::cerr
      << "Cannot create "
      << outpath
      << std::endl;

    return;
  }

  out->cd();

  TNamed(
    "note",
    "Trees are already producer-selected; original cut scan is preserved. "
    "Additional PID/angle histograms only partition candidates further.").Write();

  TNamed(
    "signed_angle_definition",
    "q_bachelor>0: phi_bachelor-phi_Lambda; "
    "q_bachelor<0: phi_Lambda-phi_bachelor.").Write();

  TNamed(
    "pid_definition",
    "Sigma/Xi bachelor PID uses pion ALEPH Nsigma; Omega uses kaon ALEPH Nsigma. "
    "Pion width uses universal width; kaon uses the momentum-dependent width correction.").Write();

  TParameter<double>("beam_x_cm", beamX).Write();
  TParameter<double>("beam_y_cm", beamY).Write();

  TParameter<double>("dedx_scale", dedxScale).Write();
  TParameter<double>("pid_loose_nsigma", pidLooseNSigma).Write();
  TParameter<double>("pid_normal_nsigma", pidNormalNSigma).Write();
  TParameter<double>("pid_tight_nsigma", pidTightNSigma).Write();

  TParameter<double>("angle_loose_at_pt0_rad", angleLooseAtZero).Write();
  TParameter<double>("angle_normal_at_pt0_rad", angleNormalAtZero).Write();
  TParameter<double>("angle_tight_at_pt0_rad", angleTightAtZero).Write();
  TParameter<double>("angle_zero_at_pt_GeV", angleZeroAtPt).Write();

  // ==========================================================================
  // Sigma(1385)
  // ==========================================================================

  TChain s(sigmaTreeName);
  const int nsfiles =
    s.Add(pattern);

  if (nsfiles > 0 &&
      s.GetEntries() > 0 &&
      has(s, "mass_sigma1385"))
  {
    const std::vector<std::string> req = {
      "is_antilambda",
      "sigma_charge",
      "bachelor_ntpc_clusters",
      "lambda_mass",
      "lambda_pair_dca",
      "lambda_dira",
      "lambda_decay_x",
      "lambda_decay_y",

      // Added angular/PID branches.
      "lambda_px",
      "lambda_py",
      "lambda_pz",
      "bachelor_px",
      "bachelor_py",
      "bachelor_pz",
      "bachelor_dedx",

      "bachelor_pt",
      "bachelor_dca_xy",
      "sigma_pt",
      "mass_sigma1385"
    };

    for (const auto& b : req)
    {
      if (!has(s, b.c_str()))
      {
        std::cerr
          << "Missing Sigma branch "
          << b
          << std::endl;
        return;
      }
    }

    Int_t anti = 0;
    Int_t charge = 0;
    UInt_t ncl = 0;

    Float_t lm = 0;
    Float_t ldca = 0;
    Float_t ldira = 0;
    Float_t lx = 0;
    Float_t ly = 0;

    Float_t lpx = 0;
    Float_t lpy = 0;
    Float_t lpz = 0;

    Float_t bpx = 0;
    Float_t bpy = 0;
    Float_t bpz = 0;
    Float_t bdedx = 0;

    Float_t bpt = 0;
    Float_t bdca = 0;
    Float_t spt = 0;
    Float_t m = 0;

    s.SetBranchAddress("is_antilambda", &anti);
    s.SetBranchAddress("sigma_charge", &charge);
    s.SetBranchAddress("bachelor_ntpc_clusters", &ncl);
    s.SetBranchAddress("lambda_mass", &lm);
    s.SetBranchAddress("lambda_pair_dca", &ldca);
    s.SetBranchAddress("lambda_dira", &ldira);
    s.SetBranchAddress("lambda_decay_x", &lx);
    s.SetBranchAddress("lambda_decay_y", &ly);

    s.SetBranchAddress("lambda_px", &lpx);
    s.SetBranchAddress("lambda_py", &lpy);
    s.SetBranchAddress("lambda_pz", &lpz);

    s.SetBranchAddress("bachelor_px", &bpx);
    s.SetBranchAddress("bachelor_py", &bpy);
    s.SetBranchAddress("bachelor_pz", &bpz);
    s.SetBranchAddress("bachelor_dedx", &bdedx);

    s.SetBranchAddress("bachelor_pt", &bpt);
    s.SetBranchAddress("bachelor_dca_xy", &bdca);
    s.SetBranchAddress("sigma_pt", &spt);
    s.SetBranchAddress("mass_sigma1385", &m);

    const std::vector<SigmaCut> cuts = {
      {
        "cut00_raw",
        "producer output",
        -1,-1,-2,-1,-1,0,-1
      },
      {
        "cut01_loose",
        "dmL<40MeV, LpairDCA<3, RL>1, bach pT>0.2, nTPC>=20",
        .040,3,-2,1,.20,20,3
      },
      {
        "cut02",
        "dmL<30MeV, LpairDCA<2, LDIRA>0.75, RL>1",
        .030,2,.75,1,.20,20,3
      },
      {
        "cut03_baseline",
        "dmL<25MeV, LpairDCA<1.5, LDIRA>0.85, RL>1.5, nTPC>=25",
        .025,1.5,.85,1.5,.20,25,2.5
      },
      {
        "cut04",
        "dmL<20MeV, LpairDCA<1, LDIRA>0.90, RL>2, nTPC>=30",
        .020,1,.90,2,.20,30,2
      },
      {
        "cut05",
        "dmL<15MeV, LpairDCA<0.8, LDIRA>0.95, RL>2, bach pT>0.25",
        .015,.8,.95,2,.25,30,1.5
      },
      {
        "cut06",
        "dmL<12MeV, LpairDCA<0.6, LDIRA>0.97, RL>2.5",
        .012,.6,.97,2.5,.25,30,1
      },
      {
        "cut07",
        "dmL<10MeV, LpairDCA<0.5, LDIRA>0.98, RL>3, bach pT>0.30",
        .010,.5,.98,3,.30,35,1
      },
      {
        "cut08",
        "dmL<8MeV, LpairDCA<0.4, LDIRA>0.99, RL>3",
        .008,.4,.99,3,.30,35,.7
      },
      {
        "cut09_tight",
        "dmL<6MeV, LpairDCA<0.3, LDIRA>0.995, RL>3",
        .006,.3,.995,3,.30,35,.5
      }
    };

    TDirectory* top =
      out->mkdir("sigma1385");

    top->cd();

    TH1I* flow =
      new TH1I(
        "h_cutflow",
        "Sigma cumulative cut flow;cut;candidates",
        cuts.size(),
        0,
        cuts.size());

    std::map<std::string, SigmaH> hs;
    std::map<std::string, PidFamily> sigmaPid;

    for (size_t i = 0; i < cuts.size(); ++i)
    {
      flow->GetXaxis()->SetBinLabel(
        i+1,
        cuts[i].name.c_str());

      auto* d =
        top->mkdir(cuts[i].name.c_str());

      d->cd();

      TNamed(
        "selection",
        cuts[i].desc.c_str()).Write();

      hs[cuts[i].name] =
        bookSigma(
          d,
          cuts[i].name);

      sigmaPid[cuts[i].name] =
        bookPidFamily(
          d,
          "sigma",
          "#Sigma(1385) #rightarrow #Lambda#pi",
          1.25,
          1.55,
          "#pi",
          "#pi");
    }

    const Long64_t n =
      maxEntries < 0
        ? s.GetEntries()
        : std::min(maxEntries, s.GetEntries());

    for (Long64_t ie = 0; ie < n; ++ie)
    {
      s.GetEntry(ie);

      const double lr =
        std::hypot(
          lx-beamX,
          ly-beamY);

      const double bachelorP =
        mag3(bpx, bpy, bpz);

      const double bachelorDedxScaled =
        bdedx * dedxScale;

      const double pionNSigma =
        targetNSigma(
          bachelorP,
          bdedx,
          kPionMass,
          false,
          dedxScale);

      const double signedDphi =
        signedLambdaBachelorDeltaPhi(
          lpx,
          lpy,
          bpx,
          bpy,
          charge);

      const double open =
        openingAngle(
          lpx,
          lpy,
          lpz,
          bpx,
          bpy,
          bpz);

      const double etaLambda =
        eta(lpx, lpy, lpz);

      const double etaBachelor =
        eta(bpx, bpy, bpz);

      const double absDeta =
        std::isfinite(etaLambda) &&
        std::isfinite(etaBachelor)
          ? std::abs(etaBachelor-etaLambda)
          : std::numeric_limits<double>::quiet_NaN();

      const double phiEtaRatio =
        std::isfinite(absDeta)
          ? std::abs(signedDphi)/(absDeta+1.e-3)
          : std::numeric_limits<double>::quiet_NaN();

      for (size_t ic = 0; ic < cuts.size(); ++ic)
      {
        const auto& cut =
          cuts[ic];

        if (!passSigma(
              cut,
              lm,
              ldca,
              ldira,
              lr,
              bpt,
              ncl,
              bdca))
        {
          continue;
        }

        flow->AddBinContent(ic+1);

        auto& h =
          hs[cut.name];

        // ---------------- original histograms ----------------
        h.mass->Fill(m);
        h.pt->Fill(spt, m);
        h.lambdaMass->Fill(lm, m);
        h.lambdaR->Fill(lr, m);
        h.lambdaDca->Fill(std::abs(ldca), m);
        h.lambdaDira->Fill(ldira, m);
        h.bachelorPt->Fill(bpt, m);
        h.bachelorDca->Fill(std::abs(bdca), m);

        if (!anti && charge > 0)
          h.lamP->Fill(m);
        else if (!anti && charge < 0)
          h.lamM->Fill(m);
        else if (anti && charge > 0)
          h.alamP->Fill(m);
        else if (anti && charge < 0)
          h.alamM->Fill(m);

        // ---------------- new PID / angular histograms ----------------
        fillPidFamily(
          sigmaPid[cut.name],
          m,
          spt,
          bachelorP,
          bachelorDedxScaled,
          pionNSigma,
          signedDphi,
          open,
          absDeta,
          phiEtaRatio,
          pidLooseNSigma,
          pidNormalNSigma,
          pidTightNSigma,
          angleLooseAtZero,
          angleNormalAtZero,
          angleTightAtZero,
          angleZeroAtPt);
      }
    }

    std::cout
      << "Sigma files "
      << nsfiles
      << ", entries "
      << s.GetEntries()
      << std::endl;
  }
  else
  {
    std::cout
      << "No usable "
      << sigmaTreeName
      << " in "
      << pattern
      << std::endl;
  }

  // ==========================================================================
  // Xi / Omega
  // ==========================================================================

  TChain c(cascadeTreeName);
  const int ncfiles =
    c.Add(pattern);

  if (ncfiles > 0 &&
      c.GetEntries() > 0 &&
      has(c, "candidate_mask"))
  {
    const std::vector<std::string> req = {
      "is_antilambda",
      "cascade_charge",
      "candidate_mask",
      "bachelor_ntpc_clusters",
      "lambda_mass",
      "lambda_pair_dca",
      "cascade_pair_dca",
      "cascade_decay_radius",
      "cascade_dira",
      "lambda_flight_from_cascade",
      "lambda_dira_from_cascade",

      // Added angular/PID branches.
      "lambda_px",
      "lambda_py",
      "lambda_pz",
      "bachelor_px",
      "bachelor_py",
      "bachelor_pz",
      "bachelor_dedx",

      "bachelor_pt",
      "bachelor_dca_xy",
      "cascade_pt",
      "mass_xi",
      "mass_omega"
    };

    for (const auto& b : req)
    {
      if (!has(c, b.c_str()))
      {
        std::cerr
          << "Missing cascade branch "
          << b
          << std::endl;
        return;
      }
    }

    Int_t anti = 0;
    Int_t charge = 0;

    UInt_t mask = 0;
    UInt_t ncl = 0;

    Float_t lm = 0;
    Float_t ldca = 0;
    Float_t cdca = 0;
    Float_t cr = 0;
    Float_t cdira = 0;
    Float_t lf = 0;
    Float_t ldira = 0;

    Float_t lpx = 0;
    Float_t lpy = 0;
    Float_t lpz = 0;

    Float_t bpx = 0;
    Float_t bpy = 0;
    Float_t bpz = 0;
    Float_t bdedx = 0;

    Float_t bpt = 0;
    Float_t bdca = 0;
    Float_t cpt = 0;
    Float_t mxi = 0;
    Float_t mom = 0;

    c.SetBranchAddress("is_antilambda", &anti);
    c.SetBranchAddress("cascade_charge", &charge);
    c.SetBranchAddress("candidate_mask", &mask);
    c.SetBranchAddress("bachelor_ntpc_clusters", &ncl);
    c.SetBranchAddress("lambda_mass", &lm);
    c.SetBranchAddress("lambda_pair_dca", &ldca);
    c.SetBranchAddress("cascade_pair_dca", &cdca);
    c.SetBranchAddress("cascade_decay_radius", &cr);
    c.SetBranchAddress("cascade_dira", &cdira);
    c.SetBranchAddress("lambda_flight_from_cascade", &lf);
    c.SetBranchAddress("lambda_dira_from_cascade", &ldira);

    c.SetBranchAddress("lambda_px", &lpx);
    c.SetBranchAddress("lambda_py", &lpy);
    c.SetBranchAddress("lambda_pz", &lpz);

    c.SetBranchAddress("bachelor_px", &bpx);
    c.SetBranchAddress("bachelor_py", &bpy);
    c.SetBranchAddress("bachelor_pz", &bpz);
    c.SetBranchAddress("bachelor_dedx", &bdedx);

    c.SetBranchAddress("bachelor_pt", &bpt);
    c.SetBranchAddress("bachelor_dca_xy", &bdca);
    c.SetBranchAddress("cascade_pt", &cpt);
    c.SetBranchAddress("mass_xi", &mxi);
    c.SetBranchAddress("mass_omega", &mom);

    const std::vector<CascadeCut> cuts = {
      {
        "cut00_raw",
        "producer output",
        -1,-1,-1,0,-1,-1,-2,-1,-2,-1
      },
      {
        "cut01_producer_like",
        "dmL<40MeV, LDCA<3, bach pT>0.2, nTPC>=20",
        .040,3,.20,20,2,.5,.75,.2,.75,-1
      },
      {
        "cut02_loose",
        "dmL<30MeV, cascadeDCA<1.5, R>0.7, DIRA>0.80",
        .030,2.5,.20,20,1.5,.7,.80,.3,.80,-1
      },
      {
        "cut03_baseline",
        "dmL<25MeV, cascadeDCA<1, R>1, DIRA>0.85, Llambda>0.5",
        .025,2,.20,25,1,1,.85,.5,.85,.02
      },
      {
        "cut04",
        "dmL<20MeV, cascadeDCA<0.8, DIRA>0.90, nTPC>=30",
        .020,1.5,.20,30,.8,1,.90,.7,.90,.03
      },
      {
        "cut05",
        "LDCA<1, cascadeDCA<0.6, R>1.5, DIRA>0.93, bach pT>0.25",
        .020,1,.25,30,.6,1.5,.93,1,.93,.05
      },
      {
        "cut06",
        "dmL<15MeV, cascadeDCA<0.5, DIRA>0.95",
        .015,.8,.25,30,.5,1.5,.95,1,.95,.05
      },
      {
        "cut07",
        "cascadeDCA<0.4, R>2, DIRA>0.97, Llambda>1.5, nTPC>=35",
        .015,.6,.30,35,.4,2,.97,1.5,.97,.08
      },
      {
        "cut08",
        "dmL<12MeV, cascadeDCA<0.3, DIRA>0.98, Llambda>2",
        .012,.5,.30,35,.3,2,.98,2,.98,.10
      },
      {
        "cut09_tight",
        "dmL<10MeV, cascadeDCA<0.2, R>2.5, DIRA>0.99",
        .010,.4,.30,35,.2,2.5,.99,2,.99,.10
      }
    };

    TDirectory* top =
      out->mkdir("cascade");

    top->cd();

    TH1I* fxi =
      new TH1I(
        "h_cutflow_xi",
        "Xi cumulative cut flow;cut;candidates",
        cuts.size(),
        0,
        cuts.size());

    TH1I* fom =
      new TH1I(
        "h_cutflow_omega",
        "Omega cumulative cut flow;cut;candidates",
        cuts.size(),
        0,
        cuts.size());

    std::map<std::string, CascadeH> hc;
    std::map<std::string, PidFamily> xiPid;
    std::map<std::string, PidFamily> omegaPid;

    for (size_t i = 0; i < cuts.size(); ++i)
    {
      fxi->GetXaxis()->SetBinLabel(
        i+1,
        cuts[i].name.c_str());

      fom->GetXaxis()->SetBinLabel(
        i+1,
        cuts[i].name.c_str());

      auto* d =
        top->mkdir(cuts[i].name.c_str());

      d->cd();

      TNamed(
        "selection",
        cuts[i].desc.c_str()).Write();

      hc[cuts[i].name] =
        bookCascade(d);

      // Xi bachelor is pion.
      xiPid[cuts[i].name] =
        bookPidFamily(
          d,
          "xi",
          "#Xi #rightarrow #Lambda#pi",
          1.25,
          1.40,
          "#pi",
          "#pi");

      // Omega bachelor is kaon.
      omegaPid[cuts[i].name] =
        bookPidFamily(
          d,
          "omega",
          "#Omega #rightarrow #Lambda K",
          1.60,
          1.75,
          "K",
          "K");
    }

    const Long64_t n =
      maxEntries < 0
        ? c.GetEntries()
        : std::min(maxEntries, c.GetEntries());

    for (Long64_t ie = 0; ie < n; ++ie)
    {
      c.GetEntry(ie);

      const bool xi =
        mask & kXi;

      const bool om =
        mask & kOmega;

      const double bachelorP =
        mag3(bpx, bpy, bpz);

      const double bachelorDedxScaled =
        bdedx * dedxScale;

      const double pionNSigma =
        targetNSigma(
          bachelorP,
          bdedx,
          kPionMass,
          false,
          dedxScale);

      const double kaonNSigma =
        targetNSigma(
          bachelorP,
          bdedx,
          kKaonMass,
          true,
          dedxScale);

      const double signedDphi =
        signedLambdaBachelorDeltaPhi(
          lpx,
          lpy,
          bpx,
          bpy,
          charge);

      const double open =
        openingAngle(
          lpx,
          lpy,
          lpz,
          bpx,
          bpy,
          bpz);

      const double etaLambda =
        eta(lpx, lpy, lpz);

      const double etaBachelor =
        eta(bpx, bpy, bpz);

      const double absDeta =
        std::isfinite(etaLambda) &&
        std::isfinite(etaBachelor)
          ? std::abs(etaBachelor-etaLambda)
          : std::numeric_limits<double>::quiet_NaN();

      const double phiEtaRatio =
        std::isfinite(absDeta)
          ? std::abs(signedDphi)/(absDeta+1.e-3)
          : std::numeric_limits<double>::quiet_NaN();

      for (size_t ic = 0; ic < cuts.size(); ++ic)
      {
        const auto& cut =
          cuts[ic];

        if (!passCascade(
              cut,
              lm,
              ldca,
              bpt,
              ncl,
              cdca,
              cr,
              cdira,
              lf,
              ldira,
              bdca))
        {
          continue;
        }

        auto& h =
          hc[cut.name];

        if (xi)
        {
          fxi->AddBinContent(ic+1);

          // ---------------- original Xi histograms ----------------
          h.xi->Fill(mxi);
          h.xiPt->Fill(cpt, mxi);
          h.xiDca->Fill(std::abs(cdca), mxi);
          h.xiR->Fill(cr, mxi);
          h.xiDira->Fill(cdira, mxi);
          h.xiLFlight->Fill(lf, mxi);
          h.xiLDira->Fill(ldira, mxi);
          h.xiBDca->Fill(std::abs(bdca), mxi);
          h.xiLmass->Fill(lm, mxi);

          if (!anti && charge < 0)
            h.xim->Fill(mxi);
          else if (anti && charge > 0)
            h.xip->Fill(mxi);

          // ---------------- new Xi PID / angular histograms ----------------
          fillPidFamily(
            xiPid[cut.name],
            mxi,
            cpt,
            bachelorP,
            bachelorDedxScaled,
            pionNSigma,
            signedDphi,
            open,
            absDeta,
            phiEtaRatio,
            pidLooseNSigma,
            pidNormalNSigma,
            pidTightNSigma,
            angleLooseAtZero,
            angleNormalAtZero,
            angleTightAtZero,
            angleZeroAtPt);
        }

        if (om)
        {
          fom->AddBinContent(ic+1);

          // ---------------- original Omega histograms ----------------
          h.om->Fill(mom);
          h.omPt->Fill(cpt, mom);
          h.omDca->Fill(std::abs(cdca), mom);
          h.omR->Fill(cr, mom);
          h.omDira->Fill(cdira, mom);
          h.omLFlight->Fill(lf, mom);
          h.omLDira->Fill(ldira, mom);
          h.omBDca->Fill(std::abs(bdca), mom);
          h.omLmass->Fill(lm, mom);

          if (!anti && charge < 0)
            h.omm->Fill(mom);
          else if (anti && charge > 0)
            h.omp->Fill(mom);

          // ---------------- new Omega PID / angular histograms ----------------
          fillPidFamily(
            omegaPid[cut.name],
            mom,
            cpt,
            bachelorP,
            bachelorDedxScaled,
            kaonNSigma,
            signedDphi,
            open,
            absDeta,
            phiEtaRatio,
            pidLooseNSigma,
            pidNormalNSigma,
            pidTightNSigma,
            angleLooseAtZero,
            angleNormalAtZero,
            angleTightAtZero,
            angleZeroAtPt);
        }
      }
    }

    std::cout
      << "Cascade files "
      << ncfiles
      << ", entries "
      << c.GetEntries()
      << std::endl;
  }
  else
  {
    std::cout
      << "No usable "
      << cascadeTreeName
      << " in "
      << pattern
      << std::endl;
  }

  out->Write();
  out->Close();

  std::cout
    << "Wrote "
    << outpath
    << std::endl;
}
