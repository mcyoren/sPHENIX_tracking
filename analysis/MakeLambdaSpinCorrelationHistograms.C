// MakeLambdaSpinCorrelationHistograms.C
//
// Pair-level Lambda/Lambda-bar/K0S/Phi correlation histogram producer for the
// TpcV0CandidateTree pairTree + likeSignPairTree used by MakeK0sPairHistograms.C.
//
// This version is designed for the EXISTING tree format produced by the
// attached TpcV0CandidateTree: run/evt, track_id1/track_id2,
// crossing_valid1/crossing_valid2, crossing1/crossing2, candidate_mask, V0
// secondary-PCA momenta, and prompt primary-PCA momenta are all already stored.
// No producer rerun is required.
//
// The design follows the analysis flow of arXiv:2506.05499 while keeping the
// detector/reconstruction selection sPHENIX-specific:
//   * build Lambda, anti-Lambda, K0S and Phi candidates from existing rows;
//   * combine two candidates from the SAME evt frame using FOUR DISTINCT tracks;
//   * fill Lambda-Lambda, Lambda-antiLambda, antiLambda-antiLambda, K0S-K0S,
//     K0S-Lambda, K0S-antiLambda, Lambda-Phi, antiLambda-Phi and K0S-Phi;
//   * save 2D double-mass and mass-mass-cos(theta*) histograms;
//   * boost the (anti)proton into each parent rest frame and form cos(theta*);
//   * split short/long range in Delta y / Delta phi;
//   * build US-US, US-LS and LS-LS background controls;
//   * build mixed-event acceptance reference matched in relative pT/phi/y;
//   * split EVERY pair histogram by crossing consistency:
//       - all pairs;
//       - both candidates have internally matching daughter crossings;
//       - both candidates internally match AND resolve to the same crossing;
//       - both candidates internally match BUT resolve to different crossings;
//       - internal daughter-crossing mismatch;
//       - invalid/missing crossing decision;
//   * treat crossing==0 as a legitimate selected crossing when crossing_valid=1;
//   * write a compact pairSummary tree so later fitting/sideband/systematic
//     choices do not require rerunning the large input trees.
//
// EXISTING TREE REQUIREMENTS
// --------------------------
// Defaults match TpcV0CandidateTree exactly:
//     evt             : Int_t frame/event sequence
//     run             : Int_t run number
//     track_id1/2     : Int_t daughter-track identities (used only within evt)
//     crossing_valid1/2, crossing1/2
//     mass_Kshort, mass_Lambda, mass_AntiLambda, mass_Phi
//     secondary-PCA momenta px1...pz2 for V0s
//     primary-PCA momenta primary_px1...primary_pz2 for Phi
//
// Entries are grouped by (input-file, run, evt). All combinations are made
// only inside that group, so different crossings inside the same evt frame are
// intentionally retained and are separated by the crossing-category folders.
//
// Mixed-event production is optional and OFF by default to keep processing and
// output size small.
//
// The final STAR-style 2D-Gaussian signal ellipse is intentionally NOT frozen
// here. The output contains the full (m1,m2,cosTheta*) information, so the
// signal mean/sigma can be fitted from the actual dataset and then used in a
// downstream extraction without hard-coding a guessed mass window.
//
// Example:
// root -l -b -q \
// 'MakeLambdaSpinCorrelationHistograms.C("/path","HITS*_V0.root","qa","lambda_spin.root")'

#include <TBranch.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace LambdaSpinCorrelation
{
constexpr double kPionMass = 0.13957039;
constexpr double kProtonMass = 0.93827208816;
constexpr double kKaonMass = 0.493677;
constexpr double kLambdaMass = 1.115683;
constexpr double kK0SMass = 0.497611;
constexpr double kTwoPi = 2.0 * TMath::Pi();

// PDG weak-decay parameters used in Eq. (1) of arXiv:2506.05499.
constexpr double kAlphaLambda = +0.747;
constexpr double kAlphaAntiLambda = -0.757;

// Must match TpcV0CandidateTree::CandidateMask used by the existing macro.
constexpr UInt_t kCandidateKShort = 1U << 0;
constexpr UInt_t kCandidateLambda = 1U << 1;
constexpr UInt_t kCandidateAntiLambda = 1U << 2;
constexpr UInt_t kCandidatePhi = 1U << 3;

enum class Species : int
{
  Lambda = 0,
  AntiLambda = 1,
  K0S = 2,
  Phi = 3
};

enum class PairType : int
{
  LambdaLambda = 0,
  LambdaAntiLambda = 1,
  AntiLambdaAntiLambda = 2,
  K0SK0S = 3,
  K0SLambda = 4,
  K0SAntiLambda = 5,
  LambdaPhi = 6,
  AntiLambdaPhi = 7,
  K0SPhi = 8
};

enum class PairSource : int
{
  USUS = 0,
  USLS = 1,
  LSLS = 2,
  MixedEvent = 3
};

enum class CrossingClass : int
{
  All = 0,
  BothCandidatesDaughtersSame = 1,
  SameResolvedCrossing = 2,
  DifferentResolvedCrossing = 3,
  InternalMismatch = 4,
  InvalidCrossing = 5
};

struct EventKey
{
  int fileIndex{0};
  int run{0};
  int event{0};

  bool operator<(const EventKey& other) const
  {
    if (fileIndex != other.fileIndex) return fileIndex < other.fileIndex;
    if (run != other.run) return run < other.run;
    return event < other.event;
  }
};

struct InputBranches
{
  Int_t event{0};
  Int_t run{0};
  Int_t trackId1{0};
  Int_t trackId2{0};
  Int_t crossingValid1{0};
  Int_t crossingValid2{0};

  Float_t massK0S{0};
  Float_t massLambda{0};
  Float_t massAntiLambda{0};
  Float_t massPhi{0};
  Float_t massP1Pi2{0};
  Float_t massPi1P2{0};
  UInt_t candidateMask{0};

  Float_t v0Pt{0};
  Float_t v0Px{0};
  Float_t v0Py{0};
  Float_t v0Pz{0};

  Float_t primaryPx1{0};
  Float_t primaryPy1{0};
  Float_t primaryPz1{0};
  Float_t primaryPx2{0};
  Float_t primaryPy2{0};
  Float_t primaryPz2{0};
  Float_t primaryPairPx{0};
  Float_t primaryPairPy{0};
  Float_t primaryPairPz{0};
  Int_t primaryPcaValid{0};
  Float_t promptPairDCA{0};
  Int_t promptPcaValid{0};

  Float_t pcaX{0};
  Float_t pcaY{0};
  Float_t pcaZ{0};
  Float_t pca1Z{0};
  Float_t pca2Z{0};

  Float_t px1{0};
  Float_t py1{0};
  Float_t pz1{0};
  Float_t px2{0};
  Float_t py2{0};
  Float_t pz2{0};

  Float_t alpha{0};
  Float_t qT{0};
  Float_t pairDCA{0};
  Float_t Lproj{0};
  Float_t dcaXY1{0};
  Float_t dcaZ1{0};
  Float_t dcaXY2{0};
  Float_t dcaZ2{0};
  Float_t charge1{0};
  Float_t charge2{0};
  Float_t quality1{0};
  Float_t quality2{0};
  Float_t dedx1{0};
  Float_t dedx2{0};
  Short_t npoints1{0};
  Short_t npoints2{0};
  Short_t crossing1{0};
  Short_t crossing2{0};
};

struct Candidate
{
  EventKey key;
  Long64_t treeEntry{-1};
  Int_t trackId1{0};
  Int_t trackId2{0};
  Species species{Species::Lambda};
  bool likeSign{false};

  double mass{0};
  double pt{0};
  double y{0};
  double phi{0};
  TLorentzVector parent;
  TLorentzVector analyzerLab;   // p / anti-p for Lambda; deterministic pion for K0S
  TLorentzVector otherDaughter;
  TVector3 analyzerRest;

  double daughterPt1{0};
  double daughterPt2{0};
  double analyzerPt{0};
  double analyzerP{0};
  double analyzerPStar{0};
  double analyzerEta{0};
  double analyzerDedx{0};
  double pionDedx{0};

  double pairDCA{0};
  double decayR{0};
  double dira{0};
  double Lproj{0};
  double alpha{0};
  double qT{0};
  double pcaZ{0};
  double deltaPcaZ{0};
  double qualityMax{0};
  int npointsMin{0};

  bool crossingValid1{false};
  bool crossingValid2{false};
  short crossing1{0};
  short crossing2{0};
};

struct EventRecord
{
  EventKey key;
  std::vector<Candidate> usLambda;
  std::vector<Candidate> usAntiLambda;
  std::vector<Candidate> usK0S;
  std::vector<Candidate> usPhi;
  std::vector<Candidate> lsLambda;
  std::vector<Candidate> lsAntiLambda;
  std::vector<Candidate> lsK0S;
  std::vector<Candidate> lsPhi;
};

struct CandidateHists
{
  TH1F* mass{nullptr};
  TH2F* massVsPt{nullptr};
  TH2F* massVsY{nullptr};
  TH2F* massVsPairDCA{nullptr};
  TH2F* massVsDecayR{nullptr};
  TH2F* massVsDira{nullptr};
  TH2F* crossing1VsCrossing2{nullptr};
  TH2F* massVsCrossingDelta{nullptr};
  TH1I* crossingValidity{nullptr};
  TH1I* crossingStatus{nullptr};
  TH2F* analyzerPtVsParentPt{nullptr};
  TH1F* analyzerPStar{nullptr};
  TH2F* analyzerDedxVsP{nullptr};
  TH2F* pionDedxVsP{nullptr};
};

struct PairHists
{
  TH1F* cosTheta{nullptr};
  TH1F* mass1{nullptr};
  TH1F* mass2{nullptr};
  TH2F* massMass{nullptr};
  TH3F* massMassCos{nullptr};
  TH1F* deltaY{nullptr};
  TH1F* deltaPhi{nullptr};
  TH1F* deltaR{nullptr};
  TH2F* deltaYDeltaPhi{nullptr};
  TH2F* pt1Pt2{nullptr};
  TH2F* y1Y2{nullptr};
  TH2F* analyzerPt1Pt2{nullptr};
  TH2F* analyzerP1P2{nullptr};
  TH2F* analyzerEta1Eta2{nullptr};
  TH2F* analyzerDedx1Dedx2{nullptr};
  TH2F* analyzerPStar1PStar2{nullptr};
  TH2F* cosVsDeltaR{nullptr};
  TH2F* cosVsDeltaY{nullptr};
  TH2F* cosVsDeltaPhi{nullptr};
  TH2F* cosVsAvgPt{nullptr};
  TH2F* cosVsAnalyzerPtMin{nullptr};
  TH2F* resolvedCrossing1Vs2{nullptr};
  TH1I* crossingCategory{nullptr};
};

struct PairSummaryRow
{
  Int_t source{0};
  Int_t pairType{0};
  Int_t crossingClass{0};
  Int_t fileIndex1{0};
  Int_t fileIndex2{0};
  Int_t run1{0};
  Int_t run2{0};
  Int_t event1{0};
  Int_t event2{0};
  Long64_t entry1{-1};
  Long64_t entry2{-1};
  Int_t trk11{0}, trk12{0}, trk21{0}, trk22{0};
  Float_t mass1{0}, mass2{0};
  Float_t pt1{0}, pt2{0};
  Float_t y1{0}, y2{0};
  Float_t phi1{0}, phi2{0};
  Float_t deltaY{0}, deltaPhi{0}, deltaR{0};
  Float_t cosThetaStar{0};
  Float_t analyzerPt1{0}, analyzerPt2{0};
  Float_t analyzerPStar1{0}, analyzerPStar2{0};
  Int_t crossingValid11{0}, crossingValid12{0}, crossingValid21{0}, crossingValid22{0};
  Short_t crossing11{0}, crossing12{0}, crossing21{0}, crossing22{0};
  Int_t candidate1DaughtersSame{0};
  Int_t candidate2DaughtersSame{0};
  Int_t bothCandidatesDaughtersSame{0};
  Int_t pairSameCrossing{0};
};

inline double wrapPhi(double x)
{
  while (x >= TMath::Pi()) x -= kTwoPi;
  while (x < -TMath::Pi()) x += kTwoPi;
  return x;
}

inline const char* speciesName(Species s)
{
  switch (s)
  {
    case Species::Lambda: return "Lambda";
    case Species::AntiLambda: return "AntiLambda";
    case Species::K0S: return "K0S";
    case Species::Phi: return "Phi";
  }
  return "Unknown";
}

inline const char* pairTypeName(PairType t)
{
  switch (t)
  {
    case PairType::LambdaLambda: return "LambdaLambda";
    case PairType::LambdaAntiLambda: return "LambdaAntiLambda";
    case PairType::AntiLambdaAntiLambda: return "AntiLambdaAntiLambda";
    case PairType::K0SK0S: return "K0SK0S";
    case PairType::K0SLambda: return "K0SLambda";
    case PairType::K0SAntiLambda: return "K0SAntiLambda";
    case PairType::LambdaPhi: return "LambdaPhi";
    case PairType::AntiLambdaPhi: return "AntiLambdaPhi";
    case PairType::K0SPhi: return "K0SPhi";
  }
  return "Unknown";
}

inline const char* sourceName(PairSource s)
{
  switch (s)
  {
    case PairSource::USUS: return "USUS";
    case PairSource::USLS: return "USLS";
    case PairSource::LSLS: return "LSLS";
    case PairSource::MixedEvent: return "MixedEvent";
  }
  return "Unknown";
}

inline const char* crossingClassName(CrossingClass c)
{
  switch (c)
  {
    case CrossingClass::All: return "all";
    case CrossingClass::BothCandidatesDaughtersSame: return "bothCandidatesDaughtersSame";
    case CrossingClass::SameResolvedCrossing: return "sameCrossing";
    case CrossingClass::DifferentResolvedCrossing: return "differentCrossing";
    case CrossingClass::InternalMismatch: return "internalMismatch";
    case CrossingClass::InvalidCrossing: return "invalidCrossing";
  }
  return "unknown";
}

inline std::vector<Candidate>& candidates(EventRecord& e, Species s, bool like)
{
  if (!like)
  {
    if (s == Species::Lambda) return e.usLambda;
    if (s == Species::AntiLambda) return e.usAntiLambda;
    if (s == Species::K0S) return e.usK0S;
    return e.usPhi;
  }
  if (s == Species::Lambda) return e.lsLambda;
  if (s == Species::AntiLambda) return e.lsAntiLambda;
  if (s == Species::K0S) return e.lsK0S;
  return e.lsPhi;
}

inline const std::vector<Candidate>& candidates(const EventRecord& e, Species s, bool like)
{
  if (!like)
  {
    if (s == Species::Lambda) return e.usLambda;
    if (s == Species::AntiLambda) return e.usAntiLambda;
    if (s == Species::K0S) return e.usK0S;
    return e.usPhi;
  }
  if (s == Species::Lambda) return e.lsLambda;
  if (s == Species::AntiLambda) return e.lsAntiLambda;
  if (s == Species::K0S) return e.lsK0S;
  return e.lsPhi;
}

inline bool fourDistinctTracks(const Candidate& a, const Candidate& b)
{
  return a.trackId1 != b.trackId1 &&
         a.trackId1 != b.trackId2 &&
         a.trackId2 != b.trackId1 &&
         a.trackId2 != b.trackId2;
}

inline bool candidateHasValidCrossing(const Candidate& c)
{
  return c.crossingValid1 && c.crossingValid2;
}

inline bool candidateDaughtersSameCrossing(const Candidate& c)
{
  return candidateHasValidCrossing(c) && c.crossing1 == c.crossing2;
}

inline CrossingClass exclusiveCrossingClass(const Candidate& a, const Candidate& b)
{
  if (!candidateHasValidCrossing(a) || !candidateHasValidCrossing(b))
    return CrossingClass::InvalidCrossing;

  if (!candidateDaughtersSameCrossing(a) || !candidateDaughtersSameCrossing(b))
    return CrossingClass::InternalMismatch;

  return a.crossing1 == b.crossing1
    ? CrossingClass::SameResolvedCrossing
    : CrossingClass::DifferentResolvedCrossing;
}

inline bool makePairType(Species a, Species b, PairType& out)
{
  if (a == Species::Lambda && b == Species::Lambda)
  {
    out = PairType::LambdaLambda; return true;
  }
  if (a == Species::AntiLambda && b == Species::AntiLambda)
  {
    out = PairType::AntiLambdaAntiLambda; return true;
  }
  if ((a == Species::Lambda && b == Species::AntiLambda) ||
      (a == Species::AntiLambda && b == Species::Lambda))
  {
    out = PairType::LambdaAntiLambda; return true;
  }
  if (a == Species::K0S && b == Species::K0S)
  {
    out = PairType::K0SK0S; return true;
  }
  if ((a == Species::K0S && b == Species::Lambda) ||
      (a == Species::Lambda && b == Species::K0S))
  {
    out = PairType::K0SLambda; return true;
  }
  if ((a == Species::K0S && b == Species::AntiLambda) ||
      (a == Species::AntiLambda && b == Species::K0S))
  {
    out = PairType::K0SAntiLambda; return true;
  }
  if ((a == Species::Lambda && b == Species::Phi) ||
      (a == Species::Phi && b == Species::Lambda))
  {
    out = PairType::LambdaPhi; return true;
  }
  if ((a == Species::AntiLambda && b == Species::Phi) ||
      (a == Species::Phi && b == Species::AntiLambda))
  {
    out = PairType::AntiLambdaPhi; return true;
  }
  if ((a == Species::K0S && b == Species::Phi) ||
      (a == Species::Phi && b == Species::K0S))
  {
    out = PairType::K0SPhi; return true;
  }
  return false;
}

inline std::pair<Species,Species> pairTypeSpecies(PairType t)
{
  switch (t)
  {
    case PairType::LambdaLambda: return {Species::Lambda,Species::Lambda};
    case PairType::LambdaAntiLambda: return {Species::Lambda,Species::AntiLambda};
    case PairType::AntiLambdaAntiLambda: return {Species::AntiLambda,Species::AntiLambda};
    case PairType::K0SK0S: return {Species::K0S,Species::K0S};
    case PairType::K0SLambda: return {Species::K0S,Species::Lambda};
    case PairType::K0SAntiLambda: return {Species::K0S,Species::AntiLambda};
    case PairType::LambdaPhi: return {Species::Lambda,Species::Phi};
    case PairType::AntiLambdaPhi: return {Species::AntiLambda,Species::Phi};
    case PairType::K0SPhi: return {Species::K0S,Species::Phi};
  }
  return {Species::Lambda,Species::Lambda};
}

inline std::pair<double,double> massRange(Species s)
{
  if (s == Species::K0S) return {0.40,0.60};
  if (s == Species::Phi) return {0.95,1.10};
  return {1.05,1.25};
}

inline void canonicalOrder(const Candidate*& a, const Candidate*& b)
{
  PairType type;
  if (!makePairType(a->species,b->species,type)) return;
  const auto axes = pairTypeSpecies(type);
  if (a->species != axes.first || b->species != axes.second)
    std::swap(a,b);
}

inline double cosThetaStar(const Candidate& a, const Candidate& b)
{
  const double ma = a.analyzerRest.Mag();
  const double mb = b.analyzerRest.Mag();
  if (!(ma > 0.) || !(mb > 0.))
    return std::numeric_limits<double>::quiet_NaN();
  return std::clamp(a.analyzerRest.Dot(b.analyzerRest)/(ma*mb), -1.0, 1.0);
}

inline bool isShortRange(const Candidate& a, const Candidate& b)
{
  return std::abs(a.y - b.y) < 0.5 &&
         std::abs(wrapPhi(a.phi - b.phi)) < TMath::Pi()/3.0;
}

inline bool isLongRange(const Candidate& a, const Candidate& b)
{
  const double dy = std::abs(a.y - b.y);
  const double dphi = std::abs(wrapPhi(a.phi - b.phi));
  return (dy >= 0.5 && dy < 2.0) ||
         (dphi >= TMath::Pi()/3.0 && dphi < TMath::Pi());
}

inline bool relativeKinematicsMatch(const Candidate& seA,
                                    const Candidate& seB,
                                    const Candidate& meA,
                                    const Candidate& meB,
                                    double maxDeltaPt,
                                    double maxDeltaPhi,
                                    double maxDeltaY)
{
  const double seDPt = std::abs(seA.pt - seB.pt);
  const double meDPt = std::abs(meA.pt - meB.pt);
  const double seDPhi = std::abs(wrapPhi(seA.phi - seB.phi));
  const double meDPhi = std::abs(wrapPhi(meA.phi - meB.phi));
  const double seDY = std::abs(seA.y - seB.y);
  const double meDY = std::abs(meA.y - meB.y);

  return std::abs(seDPt - meDPt) < maxDeltaPt &&
         std::abs(seDPhi - meDPhi) < maxDeltaPhi &&
         std::abs(seDY - meDY) < maxDeltaY;
}

inline bool branchExists(TTree* t, const char* name)
{
  return t && name && name[0] != '\0' && t->GetBranch(name) != nullptr;
}

bool bindInput(TTree* t,
               InputBranches& b,
               const char* eventBranch,
               const char* runBranch,
               const char* trackId1Branch,
               const char* trackId2Branch,
               bool requireTrackIds)
{
  if (!t) return false;

  const std::vector<const char*> required = {
    eventBranch,
    "mass_Kshort", "mass_Lambda", "mass_AntiLambda", "mass_Phi",
    "mass_P1Pi2", "mass_Pi1P2", "candidate_mask",
    "v0_pt", "v0_px", "v0_py", "v0_pz",
    "primary_px1", "primary_py1", "primary_pz1",
    "primary_px2", "primary_py2", "primary_pz2",
    "primary_pair_px", "primary_pair_py", "primary_pair_pz",
    "primary_pca_valid", "prompt_pairDCA", "prompt_pca_valid",
    "pca_x", "pca_y", "pca_z", "pca1_z", "pca2_z",
    "px1", "py1", "pz1", "px2", "py2", "pz2",
    "alpha", "qT", "pairDCA", "Lproj",
    "dca_xy1", "dca_z1", "dca_xy2", "dca_z2",
    "charge1", "charge2", "quality1", "quality2",
    "dedx_1", "dedx_2", "npoints1", "npoints2",
    "crossing_valid1", "crossing_valid2", "crossing1", "crossing2"
  };

  bool ok = true;
  for (const char* name : required)
  {
    if (!branchExists(t, name))
    {
      std::cerr << "ERROR: tree " << t->GetName()
                << " is missing required branch '" << (name ? name : "")
                << "'\n";
      ok = false;
    }
  }

  if (requireTrackIds)
  {
    for (const char* name : {trackId1Branch, trackId2Branch})
    {
      if (!branchExists(t, name))
      {
        std::cerr << "ERROR: four-distinct-track protection requires branch '"
                  << (name ? name : "") << "' in tree " << t->GetName()
                  << "\n";
        ok = false;
      }
    }
  }

  if (!ok) return false;

  t->SetBranchAddress(eventBranch, &b.event);
  if (runBranch && runBranch[0] != '\0' && branchExists(t, runBranch))
    t->SetBranchAddress(runBranch, &b.run);
  if (requireTrackIds)
  {
    t->SetBranchAddress(trackId1Branch, &b.trackId1);
    t->SetBranchAddress(trackId2Branch, &b.trackId2);
  }

  t->SetBranchAddress("mass_Kshort", &b.massK0S);
  t->SetBranchAddress("mass_Lambda", &b.massLambda);
  t->SetBranchAddress("mass_AntiLambda", &b.massAntiLambda);
  t->SetBranchAddress("mass_Phi", &b.massPhi);
  t->SetBranchAddress("mass_P1Pi2", &b.massP1Pi2);
  t->SetBranchAddress("mass_Pi1P2", &b.massPi1P2);
  t->SetBranchAddress("candidate_mask", &b.candidateMask);
  t->SetBranchAddress("v0_pt", &b.v0Pt);
  t->SetBranchAddress("v0_px", &b.v0Px);
  t->SetBranchAddress("v0_py", &b.v0Py);
  t->SetBranchAddress("v0_pz", &b.v0Pz);

  t->SetBranchAddress("primary_px1", &b.primaryPx1);
  t->SetBranchAddress("primary_py1", &b.primaryPy1);
  t->SetBranchAddress("primary_pz1", &b.primaryPz1);
  t->SetBranchAddress("primary_px2", &b.primaryPx2);
  t->SetBranchAddress("primary_py2", &b.primaryPy2);
  t->SetBranchAddress("primary_pz2", &b.primaryPz2);
  t->SetBranchAddress("primary_pair_px", &b.primaryPairPx);
  t->SetBranchAddress("primary_pair_py", &b.primaryPairPy);
  t->SetBranchAddress("primary_pair_pz", &b.primaryPairPz);
  t->SetBranchAddress("primary_pca_valid", &b.primaryPcaValid);
  t->SetBranchAddress("prompt_pairDCA", &b.promptPairDCA);
  t->SetBranchAddress("prompt_pca_valid", &b.promptPcaValid);

  t->SetBranchAddress("pca_x", &b.pcaX);
  t->SetBranchAddress("pca_y", &b.pcaY);
  t->SetBranchAddress("pca_z", &b.pcaZ);
  t->SetBranchAddress("pca1_z", &b.pca1Z);
  t->SetBranchAddress("pca2_z", &b.pca2Z);
  t->SetBranchAddress("px1", &b.px1);
  t->SetBranchAddress("py1", &b.py1);
  t->SetBranchAddress("pz1", &b.pz1);
  t->SetBranchAddress("px2", &b.px2);
  t->SetBranchAddress("py2", &b.py2);
  t->SetBranchAddress("pz2", &b.pz2);
  t->SetBranchAddress("alpha", &b.alpha);
  t->SetBranchAddress("qT", &b.qT);
  t->SetBranchAddress("pairDCA", &b.pairDCA);
  t->SetBranchAddress("Lproj", &b.Lproj);
  t->SetBranchAddress("dca_xy1", &b.dcaXY1);
  t->SetBranchAddress("dca_z1", &b.dcaZ1);
  t->SetBranchAddress("dca_xy2", &b.dcaXY2);
  t->SetBranchAddress("dca_z2", &b.dcaZ2);
  t->SetBranchAddress("charge1", &b.charge1);
  t->SetBranchAddress("charge2", &b.charge2);
  t->SetBranchAddress("quality1", &b.quality1);
  t->SetBranchAddress("quality2", &b.quality2);
  t->SetBranchAddress("dedx_1", &b.dedx1);
  t->SetBranchAddress("dedx_2", &b.dedx2);
  t->SetBranchAddress("npoints1", &b.npoints1);
  t->SetBranchAddress("npoints2", &b.npoints2);
  t->SetBranchAddress("crossing_valid1", &b.crossingValid1);
  t->SetBranchAddress("crossing_valid2", &b.crossingValid2);
  t->SetBranchAddress("crossing1", &b.crossing1);
  t->SetBranchAddress("crossing2", &b.crossing2);

  return true;
}

bool passCurrentBaseline(const InputBranches& b,
                         double beamX,
                         double beamY,
                         bool applySignedDeltaPhiCut)
{
  // Mirrors the spirit of cut03_baseline in MakeK0sPairHistograms.C.
  const double pt1 = std::hypot(b.px1, b.py1);
  const double pt2 = std::hypot(b.px2, b.py2);
  const double decayR = std::hypot(b.pcaX - beamX, b.pcaY - beamY);
  const double dx = b.pcaX - beamX;
  const double dy = b.pcaY - beamY;
  const double dz = b.pcaZ;
  const double flightMag = std::sqrt(dx*dx + dy*dy + dz*dz);
  const double pMag = std::sqrt(b.v0Px*b.v0Px + b.v0Py*b.v0Py + b.v0Pz*b.v0Pz);
  double dira = -2.0;
  if (flightMag > 0. && pMag > 0.)
    dira = (dx*b.v0Px + dy*b.v0Py + dz*b.v0Pz)/(flightMag*pMag);

  bool pass =
    std::abs(b.pcaZ) < 15.0 &&
    std::abs(b.pca1Z - b.pca2Z) < 0.50 &&
    std::min(pt1, pt2) > 0.20 &&
    decayR > 2.0 &&
    std::abs(b.alpha) < 0.99 &&
    std::abs(b.pairDCA) < 2.0 &&
    dira > 0.85 &&
    std::max(b.quality1, b.quality2) < 15.0 &&
    std::min<int>(b.npoints1, b.npoints2) > 30;

  if (!pass || !applySignedDeltaPhiCut || b.charge1*b.charge2 >= 0.)
    return pass;

  const double phi1 = std::atan2(b.py1, b.px1);
  const double phi2 = std::atan2(b.py2, b.px2);
  const double phiPositive = b.charge1 > 0 ? phi1 : phi2;
  const double phiNegative = b.charge1 > 0 ? phi2 : phi1;
  const double signedDphi = wrapPhi(phiPositive - phiNegative);
  const double threshold = 0.8 - 0.4 * std::min<double>(b.v0Pt, 2.0);
  return signedDphi >= threshold;
}

Candidate makeCandidate(const InputBranches& b,
                        const EventKey& key,
                        Long64_t entry,
                        Species species,
                        bool likeSign,
                        double beamX,
                        double beamY)
{
  Candidate c;
  c.key = key;
  c.treeEntry = entry;
  c.trackId1 = b.trackId1;
  c.trackId2 = b.trackId2;
  c.species = species;
  c.likeSign = likeSign;
  c.crossingValid1 = b.crossingValid1 != 0;
  c.crossingValid2 = b.crossingValid2 != 0;
  c.crossing1 = b.crossing1;
  c.crossing2 = b.crossing2;

  const bool isPhi = species == Species::Phi;

  const double px1 = isPhi ? b.primaryPx1 : b.px1;
  const double py1 = isPhi ? b.primaryPy1 : b.py1;
  const double pz1 = isPhi ? b.primaryPz1 : b.pz1;
  const double px2 = isPhi ? b.primaryPx2 : b.px2;
  const double py2 = isPhi ? b.primaryPy2 : b.py2;
  const double pz2 = isPhi ? b.primaryPz2 : b.pz2;

  c.daughterPt1 = std::hypot(px1, py1);
  c.daughterPt2 = std::hypot(px2, py2);
  c.pairDCA = isPhi && b.promptPcaValid
    ? std::abs(static_cast<double>(b.promptPairDCA))
    : std::abs(static_cast<double>(b.pairDCA));

  if (!isPhi)
  {
    c.decayR = std::hypot(b.pcaX - beamX, b.pcaY - beamY);
    c.Lproj = b.Lproj;
    c.alpha = b.alpha;
    c.qT = b.qT;
    c.pcaZ = b.pcaZ;
    c.deltaPcaZ = std::abs(b.pca1Z - b.pca2Z);
    c.qualityMax = std::max<double>(b.quality1, b.quality2);
    c.npointsMin = std::min<int>(b.npoints1, b.npoints2);

    const double dx = b.pcaX - beamX;
    const double dy = b.pcaY - beamY;
    const double dz = b.pcaZ;
    const double fmag = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double pmag = std::sqrt(b.v0Px*b.v0Px + b.v0Py*b.v0Py + b.v0Pz*b.v0Pz);
    c.dira = (fmag > 0. && pmag > 0.)
      ? (dx*b.v0Px + dy*b.v0Py + dz*b.v0Pz)/(fmag*pmag)
      : std::numeric_limits<double>::quiet_NaN();
  }
  else
  {
    c.decayR = std::numeric_limits<double>::quiet_NaN();
    c.dira = std::numeric_limits<double>::quiet_NaN();
    c.Lproj = std::numeric_limits<double>::quiet_NaN();
    c.alpha = std::numeric_limits<double>::quiet_NaN();
    c.qT = std::numeric_limits<double>::quiet_NaN();
    c.pcaZ = std::numeric_limits<double>::quiet_NaN();
    c.deltaPcaZ = std::numeric_limits<double>::quiet_NaN();
    c.qualityMax = std::max<double>(b.quality1, b.quality2);
    c.npointsMin = std::min<int>(b.npoints1, b.npoints2);
  }

  const bool track1Positive = b.charge1 > 0.;
  const bool track1Negative = b.charge1 < 0.;
  const bool track1HigherPt = c.daughterPt1 >= c.daughterPt2;

  int analyzerIndex = 1;
  double analyzerMass = kPionMass;
  double otherMass = kPionMass;

  if (species == Species::Lambda)
  {
    c.mass = likeSign
      ? (track1HigherPt ? b.massP1Pi2 : b.massPi1P2)
      : b.massLambda;
    analyzerMass = kProtonMass;
    otherMass = kPionMass;
    analyzerIndex = likeSign ? (track1HigherPt ? 1 : 2)
                             : (track1Positive ? 1 : 2);
  }
  else if (species == Species::AntiLambda)
  {
    c.mass = likeSign
      ? (track1HigherPt ? b.massP1Pi2 : b.massPi1P2)
      : b.massAntiLambda;
    analyzerMass = kProtonMass;
    otherMass = kPionMass;
    analyzerIndex = likeSign ? (track1HigherPt ? 1 : 2)
                             : (track1Negative ? 1 : 2);
  }
  else if (species == Species::K0S)
  {
    c.mass = b.massK0S;
    analyzerMass = kPionMass;
    otherMass = kPionMass;
    // Deterministic control analyzer: pi+ for US, higher-pT pion for LS.
    analyzerIndex = likeSign ? (track1HigherPt ? 1 : 2)
                             : (track1Positive ? 1 : 2);
  }
  else
  {
    c.mass = b.massPhi;
    analyzerMass = kKaonMass;
    otherMass = kKaonMass;
    // Deterministic control analyzer: K+ for US, higher-pT kaon for LS.
    analyzerIndex = likeSign ? (track1HigherPt ? 1 : 2)
                             : (track1Positive ? 1 : 2);
  }

  const bool analyzerIs1 = analyzerIndex == 1;
  const double apx = analyzerIs1 ? px1 : px2;
  const double apy = analyzerIs1 ? py1 : py2;
  const double apz = analyzerIs1 ? pz1 : pz2;
  const double opx = analyzerIs1 ? px2 : px1;
  const double opy = analyzerIs1 ? py2 : py1;
  const double opz = analyzerIs1 ? pz2 : pz1;

  c.analyzerDedx = analyzerIs1 ? b.dedx1 : b.dedx2;
  c.pionDedx = analyzerIs1 ? b.dedx2 : b.dedx1;
  c.analyzerLab.SetXYZM(apx, apy, apz, analyzerMass);
  c.otherDaughter.SetXYZM(opx, opy, opz, otherMass);
  c.parent = c.analyzerLab + c.otherDaughter;

  c.pt = c.parent.Pt();
  c.phi = c.parent.Phi();
  c.y = c.parent.Rapidity();
  c.analyzerPt = c.analyzerLab.Pt();
  c.analyzerP = c.analyzerLab.P();
  c.analyzerEta = c.analyzerLab.Eta();

  TLorentzVector boosted = c.analyzerLab;
  boosted.Boost(-c.parent.BoostVector());
  c.analyzerRest = boosted.Vect();
  c.analyzerPStar = boosted.P();
  return c;
}

CandidateHists bookCandidateHists(TDirectory* d, Species species, bool like)
{
  d->cd();
  CandidateHists h;
  const auto mr = massRange(species);
  const double mlo = mr.first;
  const double mhi = mr.second;
  const char* tag = like ? "LS" : "US";

  h.mass = new TH1F("h_mass", TString::Format("%s %s candidate mass;mass [GeV/c^{2}];candidates", tag, speciesName(species)), 240, mlo, mhi);
  h.massVsPt = new TH2F("h_mass_vs_pt", ";p_{T}^{candidate} [GeV/c];mass [GeV/c^{2}]", 100, 0, 10, 240, mlo, mhi);
  h.massVsY = new TH2F("h_mass_vs_y", ";y^{candidate};mass [GeV/c^{2}]", 100, -2.5, 2.5, 240, mlo, mhi);
  h.massVsPairDCA = new TH2F("h_mass_vs_pairDCA", ";|DCA_{pair}| [cm];mass [GeV/c^{2}]", 120, 0, 6, 240, mlo, mhi);
  h.massVsDecayR = new TH2F("h_mass_vs_decayR", ";R_{decay} [cm];mass [GeV/c^{2}]", 120, 0, 60, 240, mlo, mhi);
  h.massVsDira = new TH2F("h_mass_vs_DIRA", ";DIRA;mass [GeV/c^{2}]", 120, -0.2, 1.0, 240, mlo, mhi);

  h.crossing1VsCrossing2 = new TH2F(
      "h_crossing1_vs_crossing2",
      ";daughter crossing_{1};daughter crossing_{2}",
      600, -120, 480, 600, -120, 480);
  h.massVsCrossingDelta = new TH2F(
      "h_mass_vs_delta_crossing",
      ";crossing_{1}-crossing_{2};mass [GeV/c^{2}]",
      601, -300.5, 300.5, 240, mlo, mhi);

  h.crossingValidity = new TH1I("h_crossing_validity", ";crossing validity;candidates", 4, 0, 4);
  h.crossingValidity->GetXaxis()->SetBinLabel(1, "both valid");
  h.crossingValidity->GetXaxis()->SetBinLabel(2, "only daughter1 valid");
  h.crossingValidity->GetXaxis()->SetBinLabel(3, "only daughter2 valid");
  h.crossingValidity->GetXaxis()->SetBinLabel(4, "neither valid");

  h.crossingStatus = new TH1I("h_crossing_status", ";status;candidates", 3, 0, 3);
  h.crossingStatus->GetXaxis()->SetBinLabel(1, "valid + same");
  h.crossingStatus->GetXaxis()->SetBinLabel(2, "valid + different");
  h.crossingStatus->GetXaxis()->SetBinLabel(3, "invalid");

  h.analyzerPtVsParentPt = new TH2F("h_analyzer_pt_vs_parent_pt", ";p_{T}^{candidate};p_{T}^{analyzer}", 100, 0, 10, 100, 0, 10);
  h.analyzerPStar = new TH1F("h_analyzer_pstar", ";p^{*}_{analyzer} [GeV/c];candidates", 150, 0, 1.5);
  h.analyzerDedxVsP = new TH2F("h_analyzer_dedx_vs_p", ";p_{analyzer} [GeV/c];dE/dx", 160, 0, 8, 180, 0, 900);
  h.pionDedxVsP = new TH2F("h_other_daughter_dedx_vs_p", ";p_{other} [GeV/c];dE/dx", 160, 0, 8, 180, 0, 900);
  return h;
}

PairHists bookPairHists(TDirectory* d, PairType type)
{
  d->cd();
  PairHists h;
  const auto axes = pairTypeSpecies(type);
  const auto mr1 = massRange(axes.first);
  const auto mr2 = massRange(axes.second);

  h.cosTheta = new TH1F("h_cosThetaStar", ";cos#theta^{*};pairs", 40, -1, 1);
  h.mass1 = new TH1F("h_mass1", TString::Format(";m_{1}(%s) [GeV/c^{2}];pairs", speciesName(axes.first)), 240, mr1.first, mr1.second);
  h.mass2 = new TH1F("h_mass2", TString::Format(";m_{2}(%s) [GeV/c^{2}];pairs", speciesName(axes.second)), 240, mr2.first, mr2.second);
  h.massMass = new TH2F("h_mass1_vs_mass2",
      TString::Format(";m_{1}(%s) [GeV/c^{2}];m_{2}(%s) [GeV/c^{2}]",
                      speciesName(axes.first),speciesName(axes.second)),
      120, mr1.first, mr1.second, 120, mr2.first, mr2.second);
  h.massMassCos = new TH3F("h3_mass1_mass2_cosThetaStar",
      TString::Format(";m_{1}(%s);m_{2}(%s);cos#theta^{*}",
                      speciesName(axes.first),speciesName(axes.second)),
      80, mr1.first, mr1.second, 80, mr2.first, mr2.second, 30, -1, 1);
  h.deltaY = new TH1F("h_absDeltaY", ";|#Deltay|;pairs", 100, 0, 2.5);
  h.deltaPhi = new TH1F("h_absDeltaPhi", ";|#Delta#phi|;pairs", 100, 0, TMath::Pi());
  h.deltaR = new TH1F("h_deltaR", ";#DeltaR=#sqrt{(#Deltay)^{2}+(#Delta#phi)^{2}};pairs", 120, 0, 4.0);
  h.deltaYDeltaPhi = new TH2F("h_absDeltaY_vs_absDeltaPhi", ";|#Deltay|;|#Delta#phi|", 80, 0, 2.0, 80, 0, TMath::Pi());
  h.pt1Pt2 = new TH2F("h_pt1_vs_pt2", ";p_{T,1};p_{T,2}", 80, 0, 8, 80, 0, 8);
  h.y1Y2 = new TH2F("h_y1_vs_y2", ";y_{1};y_{2}", 80, -1.2, 1.2, 80, -1.2, 1.2);
  h.analyzerPt1Pt2 = new TH2F("h_analyzer_pt1_vs_pt2", ";p_{T,1}^{analyzer};p_{T,2}^{analyzer}", 100, 0, 6, 100, 0, 6);
  h.analyzerP1P2 = new TH2F("h_analyzer_p1_vs_p2", ";p_{1}^{analyzer};p_{2}^{analyzer}", 100, 0, 8, 100, 0, 8);
  h.analyzerEta1Eta2 = new TH2F("h_analyzer_eta1_vs_eta2", ";#eta_{1}^{analyzer};#eta_{2}^{analyzer}", 100, -2, 2, 100, -2, 2);
  h.analyzerDedx1Dedx2 = new TH2F("h_analyzer_dedx1_vs_dedx2", ";dE/dx_{1}^{analyzer};dE/dx_{2}^{analyzer}", 120, 0, 900, 120, 0, 900);
  h.analyzerPStar1PStar2 = new TH2F("h_analyzer_pstar1_vs_pstar2", ";p^{*}_{1};p^{*}_{2}", 100, 0, 1, 100, 0, 1);
  h.cosVsDeltaR = new TH2F("h_cosThetaStar_vs_deltaR", ";#DeltaR;cos#theta^{*}", 80, 0, 4.0, 40, -1, 1);
  h.cosVsDeltaY = new TH2F("h_cosThetaStar_vs_absDeltaY", ";|#Deltay|;cos#theta^{*}", 80, 0, 2.5, 40, -1, 1);
  h.cosVsDeltaPhi = new TH2F("h_cosThetaStar_vs_absDeltaPhi", ";|#Delta#phi|;cos#theta^{*}", 80, 0, TMath::Pi(), 40, -1, 1);
  h.cosVsAvgPt = new TH2F("h_cosThetaStar_vs_avg_pair_pt", ";(p_{T,1}+p_{T,2})/2;cos#theta^{*}", 80, 0, 6, 40, -1, 1);
  h.cosVsAnalyzerPtMin = new TH2F("h_cosThetaStar_vs_min_analyzer_pt", ";min p_{T}^{analyzer};cos#theta^{*}", 80, 0, 4, 40, -1, 1);
  h.resolvedCrossing1Vs2 = new TH2F("h_resolved_crossing1_vs_crossing2", ";candidate crossing 1;candidate crossing 2", 600, -120, 480, 600, -120, 480);
  h.crossingCategory = new TH1I("h_crossing_category", ";crossing category;pairs", 4, 0, 4);
  h.crossingCategory->GetXaxis()->SetBinLabel(1, "same crossing");
  h.crossingCategory->GetXaxis()->SetBinLabel(2, "different crossing");
  h.crossingCategory->GetXaxis()->SetBinLabel(3, "internal mismatch");
  h.crossingCategory->GetXaxis()->SetBinLabel(4, "invalid");
  return h;
}

void fillCandidateHists(CandidateHists& h, const Candidate& c)
{
  h.mass->Fill(c.mass);
  h.massVsPt->Fill(c.pt, c.mass);
  h.massVsY->Fill(c.y, c.mass);
  if (std::isfinite(c.pairDCA)) h.massVsPairDCA->Fill(c.pairDCA, c.mass);
  if (std::isfinite(c.decayR)) h.massVsDecayR->Fill(c.decayR, c.mass);
  if (std::isfinite(c.dira)) h.massVsDira->Fill(c.dira, c.mass);

  if (c.crossingValid1 && c.crossingValid2)
  {
    h.crossingValidity->Fill(0.5);
    h.crossing1VsCrossing2->Fill(c.crossing1, c.crossing2);
    h.massVsCrossingDelta->Fill(static_cast<int>(c.crossing1)-static_cast<int>(c.crossing2), c.mass);
    h.crossingStatus->Fill(c.crossing1 == c.crossing2 ? 0.5 : 1.5);
  }
  else
  {
    if (c.crossingValid1 && !c.crossingValid2) h.crossingValidity->Fill(1.5);
    else if (!c.crossingValid1 && c.crossingValid2) h.crossingValidity->Fill(2.5);
    else h.crossingValidity->Fill(3.5);
    h.crossingStatus->Fill(2.5);
  }

  h.analyzerPtVsParentPt->Fill(c.pt, c.analyzerPt);
  h.analyzerPStar->Fill(c.analyzerPStar);
  if (std::isfinite(c.analyzerDedx)) h.analyzerDedxVsP->Fill(c.analyzerP, c.analyzerDedx);
  if (std::isfinite(c.pionDedx)) h.pionDedxVsP->Fill(c.otherDaughter.P(), c.pionDedx);
}

void fillPairHists(PairHists& h, const Candidate& inputA, const Candidate& inputB, double weight=1.0)
{
  const Candidate* a = &inputA;
  const Candidate* b = &inputB;
  canonicalOrder(a, b);
  const double c = cosThetaStar(*a, *b);
  if (!std::isfinite(c)) return;
  const double dy = std::abs(a->y - b->y);
  const double dphi = std::abs(wrapPhi(a->phi - b->phi));
  const double dr = std::hypot(dy, dphi);

  h.cosTheta->Fill(c, weight);
  h.mass1->Fill(a->mass, weight);
  h.mass2->Fill(b->mass, weight);
  h.massMass->Fill(a->mass, b->mass, weight);
  h.massMassCos->Fill(a->mass, b->mass, c, weight);
  h.deltaY->Fill(dy, weight);
  h.deltaPhi->Fill(dphi, weight);
  h.deltaR->Fill(dr, weight);
  h.deltaYDeltaPhi->Fill(dy, dphi, weight);
  h.pt1Pt2->Fill(a->pt, b->pt, weight);
  h.y1Y2->Fill(a->y, b->y, weight);
  h.analyzerPt1Pt2->Fill(a->analyzerPt, b->analyzerPt, weight);
  h.analyzerP1P2->Fill(a->analyzerP, b->analyzerP, weight);
  h.analyzerEta1Eta2->Fill(a->analyzerEta, b->analyzerEta, weight);
  if (std::isfinite(a->analyzerDedx) && std::isfinite(b->analyzerDedx))
    h.analyzerDedx1Dedx2->Fill(a->analyzerDedx, b->analyzerDedx, weight);
  h.analyzerPStar1PStar2->Fill(a->analyzerPStar, b->analyzerPStar, weight);
  h.cosVsDeltaR->Fill(dr, c, weight);
  h.cosVsDeltaY->Fill(dy, c, weight);
  h.cosVsDeltaPhi->Fill(dphi, c, weight);
  h.cosVsAvgPt->Fill(0.5*(a->pt+b->pt), c, weight);
  h.cosVsAnalyzerPtMin->Fill(std::min(a->analyzerPt,b->analyzerPt), c, weight);

  const CrossingClass cc = exclusiveCrossingClass(*a,*b);
  if (cc == CrossingClass::SameResolvedCrossing)
  {
    h.crossingCategory->Fill(0.5,weight);
    h.resolvedCrossing1Vs2->Fill(a->crossing1,b->crossing1,weight);
  }
  else if (cc == CrossingClass::DifferentResolvedCrossing)
  {
    h.crossingCategory->Fill(1.5,weight);
    h.resolvedCrossing1Vs2->Fill(a->crossing1,b->crossing1,weight);
  }
  else if (cc == CrossingClass::InternalMismatch)
  {
    h.crossingCategory->Fill(2.5,weight);
  }
  else
  {
    h.crossingCategory->Fill(3.5,weight);
  }
}

struct HistogramStore
{
  std::map<std::string, CandidateHists> candidate;
  std::map<std::string, PairHists> pair;

  static std::string candidateKey(Species s, bool like)
  {
    return std::string(like ? "LS/" : "US/") + speciesName(s);
  }

  static std::string pairKey(PairSource source, PairType type,
                             const std::string& range, CrossingClass crossing)
  {
    return std::string(sourceName(source)) + "/" + pairTypeName(type) +
           "/" + range + "/" + crossingClassName(crossing);
  }
};

void fillPairWithCategories(HistogramStore& hs,
                            PairSource source,
                            const Candidate& a,
                            const Candidate& b,
                            double weight=1.0)
{
  PairType type;
  if (!makePairType(a.species, b.species, type)) return;

  const CrossingClass exclusive = exclusiveCrossingClass(a,b);
  const bool bothDaughtersSame =
      candidateDaughtersSameCrossing(a) &&
      candidateDaughtersSameCrossing(b);

  const bool shortRange = isShortRange(a, b);
  const bool longRange = isLongRange(a, b);

  for (const std::string& range : {std::string("all"),
                                   shortRange ? std::string("short") : std::string(),
                                   longRange ? std::string("long") : std::string()})
  {
    if (range.empty()) continue;

    fillPairHists(
        hs.pair.at(HistogramStore::pairKey(
            source,type,range,CrossingClass::All)),
        a,b,weight);

    // Crossing equality is meaningful only for candidates reconstructed in the
    // same evt frame. Mixed-event pairs therefore fill only the inclusive
    // crossing folder.
    if (source == PairSource::MixedEvent)
      continue;

    if (bothDaughtersSame)
    {
      fillPairHists(
          hs.pair.at(HistogramStore::pairKey(
              source,type,range,CrossingClass::BothCandidatesDaughtersSame)),
          a,b,weight);
    }

    fillPairHists(
        hs.pair.at(HistogramStore::pairKey(
            source,type,range,exclusive)),
        a,b,weight);
  }
}

void setupPairSummaryTree(TTree* t, PairSummaryRow& r)
{
  t->Branch("source", &r.source, "source/I");
  t->Branch("pair_type", &r.pairType, "pair_type/I");
  t->Branch("crossing_class", &r.crossingClass, "crossing_class/I");
  t->Branch("file_index1", &r.fileIndex1, "file_index1/I");
  t->Branch("file_index2", &r.fileIndex2, "file_index2/I");
  t->Branch("run1", &r.run1, "run1/I");
  t->Branch("run2", &r.run2, "run2/I");
  t->Branch("event1", &r.event1, "event1/I");
  t->Branch("event2", &r.event2, "event2/I");
  t->Branch("entry1", &r.entry1, "entry1/L");
  t->Branch("entry2", &r.entry2, "entry2/L");
  t->Branch("track11", &r.trk11, "track11/I");
  t->Branch("track12", &r.trk12, "track12/I");
  t->Branch("track21", &r.trk21, "track21/I");
  t->Branch("track22", &r.trk22, "track22/I");
  t->Branch("mass1", &r.mass1, "mass1/F");
  t->Branch("mass2", &r.mass2, "mass2/F");
  t->Branch("pt1", &r.pt1, "pt1/F");
  t->Branch("pt2", &r.pt2, "pt2/F");
  t->Branch("y1", &r.y1, "y1/F");
  t->Branch("y2", &r.y2, "y2/F");
  t->Branch("phi1", &r.phi1, "phi1/F");
  t->Branch("phi2", &r.phi2, "phi2/F");
  t->Branch("delta_y", &r.deltaY, "delta_y/F");
  t->Branch("delta_phi", &r.deltaPhi, "delta_phi/F");
  t->Branch("delta_r", &r.deltaR, "delta_r/F");
  t->Branch("cos_theta_star", &r.cosThetaStar, "cos_theta_star/F");
  t->Branch("analyzer_pt1", &r.analyzerPt1, "analyzer_pt1/F");
  t->Branch("analyzer_pt2", &r.analyzerPt2, "analyzer_pt2/F");
  t->Branch("analyzer_pstar1", &r.analyzerPStar1, "analyzer_pstar1/F");
  t->Branch("analyzer_pstar2", &r.analyzerPStar2, "analyzer_pstar2/F");

  t->Branch("crossing_valid11", &r.crossingValid11, "crossing_valid11/I");
  t->Branch("crossing_valid12", &r.crossingValid12, "crossing_valid12/I");
  t->Branch("crossing_valid21", &r.crossingValid21, "crossing_valid21/I");
  t->Branch("crossing_valid22", &r.crossingValid22, "crossing_valid22/I");
  t->Branch("crossing11", &r.crossing11, "crossing11/S");
  t->Branch("crossing12", &r.crossing12, "crossing12/S");
  t->Branch("crossing21", &r.crossing21, "crossing21/S");
  t->Branch("crossing22", &r.crossing22, "crossing22/S");
  t->Branch("candidate1_daughters_same_crossing",
            &r.candidate1DaughtersSame,
            "candidate1_daughters_same_crossing/I");
  t->Branch("candidate2_daughters_same_crossing",
            &r.candidate2DaughtersSame,
            "candidate2_daughters_same_crossing/I");
  t->Branch("both_candidates_daughters_same_crossing",
            &r.bothCandidatesDaughtersSame,
            "both_candidates_daughters_same_crossing/I");
  t->Branch("pair_same_crossing",
            &r.pairSameCrossing,
            "pair_same_crossing/I");
}

void fillPairSummary(TTree* t, PairSummaryRow& r,
                     PairSource source,
                     const Candidate& inputA,
                     const Candidate& inputB)
{
  if (!t) return;
  const Candidate* a = &inputA;
  const Candidate* b = &inputB;
  canonicalOrder(a, b);

  PairType type;
  if (!makePairType(a->species,b->species,type)) return;

  const double c = cosThetaStar(*a, *b);
  if (!std::isfinite(c)) return;
  const double dy = std::abs(a->y-b->y);
  const double dphi = std::abs(wrapPhi(a->phi-b->phi));

  const bool aSame = candidateDaughtersSameCrossing(*a);
  const bool bSame = candidateDaughtersSameCrossing(*b);
  const bool bothSame = aSame && bSame;
  const bool pairSame = bothSame && a->crossing1 == b->crossing1;

  r.source = static_cast<int>(source);
  r.pairType = static_cast<int>(type);
  r.crossingClass = static_cast<int>(exclusiveCrossingClass(*a,*b));
  r.fileIndex1 = a->key.fileIndex;
  r.fileIndex2 = b->key.fileIndex;
  r.run1 = a->key.run; r.run2 = b->key.run;
  r.event1 = a->key.event; r.event2 = b->key.event;
  r.entry1 = a->treeEntry; r.entry2 = b->treeEntry;
  r.trk11 = a->trackId1; r.trk12 = a->trackId2;
  r.trk21 = b->trackId1; r.trk22 = b->trackId2;
  r.mass1 = a->mass; r.mass2 = b->mass;
  r.pt1 = a->pt; r.pt2 = b->pt;
  r.y1 = a->y; r.y2 = b->y;
  r.phi1 = a->phi; r.phi2 = b->phi;
  r.deltaY = dy; r.deltaPhi = dphi; r.deltaR = std::hypot(dy,dphi);
  r.cosThetaStar = c;
  r.analyzerPt1 = a->analyzerPt; r.analyzerPt2 = b->analyzerPt;
  r.analyzerPStar1 = a->analyzerPStar; r.analyzerPStar2 = b->analyzerPStar;

  r.crossingValid11 = a->crossingValid1 ? 1 : 0;
  r.crossingValid12 = a->crossingValid2 ? 1 : 0;
  r.crossingValid21 = b->crossingValid1 ? 1 : 0;
  r.crossingValid22 = b->crossingValid2 ? 1 : 0;
  r.crossing11 = a->crossing1; r.crossing12 = a->crossing2;
  r.crossing21 = b->crossing1; r.crossing22 = b->crossing2;
  r.candidate1DaughtersSame = aSame ? 1 : 0;
  r.candidate2DaughtersSame = bSame ? 1 : 0;
  r.bothCandidatesDaughtersSame = bothSame ? 1 : 0;
  r.pairSameCrossing = pairSame ? 1 : 0;

  t->Fill();
}

void pairSameSpecies(const std::vector<Candidate>& v,
                     PairSource source,
                     HistogramStore& hs,
                     TTree* summary,
                     PairSummaryRow& row,
                     bool requireFourDistinct)
{
  for (std::size_t i=0;i<v.size();++i)
  {
    for (std::size_t j=i+1;j<v.size();++j)
    {
      if (requireFourDistinct && !fourDistinctTracks(v[i],v[j])) continue;
      fillPairWithCategories(hs, source, v[i], v[j]);
      fillPairSummary(summary,row,source,v[i],v[j]);
    }
  }
}

void pairTwoSpecies(const std::vector<Candidate>& a,
                    const std::vector<Candidate>& b,
                    PairSource source,
                    HistogramStore& hs,
                    TTree* summary,
                    PairSummaryRow& row,
                    bool requireFourDistinct)
{
  for (const auto& x : a)
  {
    for (const auto& y : b)
    {
      if (requireFourDistinct && !fourDistinctTracks(x,y)) continue;
      fillPairWithCategories(hs, source, x, y);
      fillPairSummary(summary,row,source,x,y);
    }
  }
}

void pairUSLS(const EventRecord& e,
              HistogramStore& hs,
              TTree* summary,
              PairSummaryRow& row,
              bool requireFourDistinct)
{
  // Both orientations are included for unlike-species combinations so US-LS
  // control samples remain symmetric with respect to which leg is the real
  // unlike-sign candidate.

  pairTwoSpecies(e.usLambda, e.lsLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usAntiLambda, e.lsAntiLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usK0S, e.lsK0S, PairSource::USLS, hs, summary, row, requireFourDistinct);

  pairTwoSpecies(e.usLambda, e.lsAntiLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsLambda, e.usAntiLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);

  pairTwoSpecies(e.usK0S, e.lsLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsK0S, e.usLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);

  pairTwoSpecies(e.usK0S, e.lsAntiLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsK0S, e.usAntiLambda, PairSource::USLS, hs, summary, row, requireFourDistinct);

  pairTwoSpecies(e.usLambda, e.lsPhi, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsLambda, e.usPhi, PairSource::USLS, hs, summary, row, requireFourDistinct);

  pairTwoSpecies(e.usAntiLambda, e.lsPhi, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsAntiLambda, e.usPhi, PairSource::USLS, hs, summary, row, requireFourDistinct);

  pairTwoSpecies(e.usK0S, e.lsPhi, PairSource::USLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsK0S, e.usPhi, PairSource::USLS, hs, summary, row, requireFourDistinct);
}

void processSameEvent(const EventRecord& e,
                      HistogramStore& hs,
                      TTree* summary,
                      PairSummaryRow& row,
                      bool requireFourDistinct)
{
  // US-US signal/control channels.
  pairSameSpecies(e.usLambda, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairSameSpecies(e.usAntiLambda, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairSameSpecies(e.usK0S, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usLambda, e.usAntiLambda, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usK0S, e.usLambda, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usK0S, e.usAntiLambda, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usLambda, e.usPhi, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usAntiLambda, e.usPhi, PairSource::USUS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.usK0S, e.usPhi, PairSource::USUS, hs, summary, row, requireFourDistinct);

  // US-LS background controls.
  pairUSLS(e, hs, summary, row, requireFourDistinct);

  // LS-LS controls.
  pairSameSpecies(e.lsLambda, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairSameSpecies(e.lsAntiLambda, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairSameSpecies(e.lsK0S, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsLambda, e.lsAntiLambda, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsK0S, e.lsLambda, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsK0S, e.lsAntiLambda, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsLambda, e.lsPhi, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsAntiLambda, e.lsPhi, PairSource::LSLS, hs, summary, row, requireFourDistinct);
  pairTwoSpecies(e.lsK0S, e.lsPhi, PairSource::LSLS, hs, summary, row, requireFourDistinct);
}

void mixedForSEPair(const Candidate& seA,
                    const Candidate& seB,
                    const std::deque<EventRecord>& pool,
                    HistogramStore& hs,
                    double maxDeltaPt,
                    double maxDeltaPhi,
                    double maxDeltaY)
{
  std::vector<std::pair<const Candidate*,const Candidate*>> matches;
  for (const auto& old : pool)
  {
    const auto& bCandidates = candidates(old, seB.species, false);
    for (const auto& mixB : bCandidates)
    {
      if (relativeKinematicsMatch(seA,seB,seA,mixB,maxDeltaPt,maxDeltaPhi,maxDeltaY))
        matches.emplace_back(&seA,&mixB);
    }

    const auto& aCandidates = candidates(old, seA.species, false);
    for (const auto& mixA : aCandidates)
    {
      if (relativeKinematicsMatch(seA,seB,mixA,seB,maxDeltaPt,maxDeltaPhi,maxDeltaY))
        matches.emplace_back(&mixA,&seB);
    }
  }

  if (matches.empty()) return;
  const double weight = 1.0/static_cast<double>(matches.size());
  for (const auto& m : matches)
    fillPairWithCategories(hs, PairSource::MixedEvent, *m.first, *m.second, weight);
}

void buildMixedEvent(const EventRecord& e,
                     const std::deque<EventRecord>& pool,
                     HistogramStore& hs,
                     bool requireFourDistinct,
                     double maxDeltaPt,
                     double maxDeltaPhi,
                     double maxDeltaY)
{
  auto same = [&](const std::vector<Candidate>& v)
  {
    for (std::size_t i=0;i<v.size();++i)
      for (std::size_t j=i+1;j<v.size();++j)
      {
        if (requireFourDistinct && !fourDistinctTracks(v[i],v[j])) continue;
        mixedForSEPair(v[i],v[j],pool,hs,maxDeltaPt,maxDeltaPhi,maxDeltaY);
      }
  };

  auto cross = [&](const std::vector<Candidate>& a,
                   const std::vector<Candidate>& b)
  {
    for (const auto& x : a)
      for (const auto& y : b)
      {
        if (requireFourDistinct && !fourDistinctTracks(x,y)) continue;
        mixedForSEPair(x,y,pool,hs,maxDeltaPt,maxDeltaPhi,maxDeltaY);
      }
  };

  same(e.usLambda);
  same(e.usAntiLambda);
  same(e.usK0S);
  cross(e.usLambda,e.usAntiLambda);
  cross(e.usK0S,e.usLambda);
  cross(e.usK0S,e.usAntiLambda);
  cross(e.usLambda,e.usPhi);
  cross(e.usAntiLambda,e.usPhi);
  cross(e.usK0S,e.usPhi);
}

void addCandidatesFromTree(TTree* t,
                           bool likeSign,
                           int fileIndex,
                           std::map<EventKey,EventRecord>& events,
                           HistogramStore& hs,
                           Long64_t& processedEntries,
                           Long64_t maxEntries,
                           const char* eventBranch,
                           const char* runBranch,
                           const char* trackId1Branch,
                           const char* trackId2Branch,
                           bool requireTrackIds,
                           bool applyCurrentSignedDeltaPhiCut,
                           bool requireProtonHigherPt,
                           double beamX,
                           double beamY)
{
  if (!t) return;

  InputBranches b;
  if (!bindInput(t,b,eventBranch,runBranch,trackId1Branch,trackId2Branch,
                 requireTrackIds))
    return;

  const Long64_t n = t->GetEntries();
  for (Long64_t i=0;i<n;++i)
  {
    if (maxEntries >= 0 && processedEntries >= maxEntries) break;

    t->GetEntry(i);
    ++processedEntries;

    const bool unlike = b.charge1*b.charge2 < 0.;
    const bool like = b.charge1*b.charge2 > 0.;
    if ((!likeSign && !unlike) || (likeSign && !like)) continue;

    EventKey key;
    key.fileIndex = fileIndex;
    key.run = (runBranch && runBranch[0] != '\0' && branchExists(t,runBranch))
      ? b.run : 0;
    key.event = b.event;

    auto& event = events[key];
    event.key = key;

    const bool v0BaselinePass =
      passCurrentBaseline(b,beamX,beamY,
                          likeSign ? false : applyCurrentSignedDeltaPhiCut);

    auto acceptKinematics = [](const Candidate& c)
    {
      return std::isfinite(c.y) &&
             std::abs(c.y) < 1.0 &&
             c.pt > 0.5 &&
             c.pt < 5.0;
    };

    auto broadMass = [](Species s, double mass)
    {
      if (!std::isfinite(mass)) return false;
      if (s == Species::K0S) return mass > 0.40 && mass < 0.60;
      if (s == Species::Phi) return mass > 0.95 && mass < 1.10;
      return mass > 1.05 && mass < 1.25;
    };

    auto protonHigherPtPass = [&](Species s)
    {
      if (!requireProtonHigherPt || likeSign ||
          (s != Species::Lambda && s != Species::AntiLambda))
        return true;

      const double pt1 = std::hypot(b.px1,b.py1);
      const double pt2 = std::hypot(b.px2,b.py2);

      if (s == Species::Lambda)
        return b.charge1 > 0 ? pt1 > pt2 : pt2 > pt1;

      return b.charge1 < 0 ? pt1 > pt2 : pt2 > pt1;
    };

    std::vector<Species> toBuild;

    // V0 species: keep the offline V0 baseline used by this macro in addition
    // to the producer candidate_mask.
    if (v0BaselinePass)
    {
      if (b.candidateMask & kCandidateLambda)
        toBuild.push_back(Species::Lambda);
      if (b.candidateMask & kCandidateAntiLambda)
        toBuild.push_back(Species::AntiLambda);
      if (b.candidateMask & kCandidateKShort)
        toBuild.push_back(Species::K0S);
    }

    // Phi is prompt and must NOT be subjected to displaced-V0 cuts. The
    // producer already evaluated its prompt selection and stored primary-PCA
    // daughter momenta. The candidate mask is therefore the authoritative
    // inclusion flag here.
    if ((b.candidateMask & kCandidatePhi) && b.primaryPcaValid)
      toBuild.push_back(Species::Phi);

    for (Species species : toBuild)
    {
      if (!protonHigherPtPass(species)) continue;

      Candidate c =
        makeCandidate(b,key,i,species,likeSign,beamX,beamY);

      if (!broadMass(species,c.mass) || !acceptKinematics(c))
        continue;

      candidates(event,species,likeSign).push_back(c);
      fillCandidateHists(
          hs.candidate.at(
              HistogramStore::candidateKey(species,likeSign)),
          c);
    }
  }
}

} // namespace LambdaSpinCorrelation

void MakeLambdaSpinCorrelationHistograms(
  const char* inputDir = ".",
  const char* filePattern = "*.root",
  const char* outputDir = "output",
  const char* outputName = "lambda_spin_correlations.root",
  const char* pairTreeName = "pairTree",
  const char* likeSignTreeName = "likeSignPairTree",
  const char* eventBranch = "evt",
  const char* runBranch = "run",
  const char* trackId1Branch = "track_id1",
  const char* trackId2Branch = "track_id2",
  const bool requireFourDistinctTracks = true,
  const bool includeLikeSignTree = true,
  const bool doMixedEvent = false,
  const int mixDepth = 5,
  const bool applyCurrentSignedDeltaPhiCut = true,
  const bool requireProtonHigherPt = false,
  const double beamX = 0.158,
  const double beamY = 0.285,
  const Long64_t maxEntries = -1,
  const bool writePairSummary = true,
  const double mixedMatchDeltaPt = 0.10,
  const double mixedMatchDeltaPhi = 0.10,
  const double mixedMatchDeltaY = 0.10)
{
  using namespace LambdaSpinCorrelation;
  TH1::AddDirectory(kTRUE);

  if (mixDepth < 1 && doMixedEvent)
  {
    std::cerr << "ERROR: mixDepth must be >=1 when doMixedEvent=true\n";
    return;
  }

  const TString inputPattern =
      TString::Format("%s/%s",inputDir,filePattern);

  TChain fileFinder(pairTreeName);
  const int nFiles = fileFinder.Add(inputPattern);
  if (nFiles <= 0)
  {
    std::cerr << "ERROR: no files matched "
              << inputPattern << "\n";
    return;
  }

  std::vector<std::string> files;
  if (auto* list = fileFinder.GetListOfFiles())
  {
    for (int i=0;i<list->GetEntries();++i)
    {
      auto* el =
        dynamic_cast<TChainElement*>(list->At(i));
      if (el) files.emplace_back(el->GetTitle());
    }
  }

  gSystem->mkdir(outputDir,kTRUE);
  const TString outputPath =
      TString::Format("%s/%s",outputDir,outputName);

  std::unique_ptr<TFile> out(
      TFile::Open(outputPath,"RECREATE"));

  if (!out || out->IsZombie())
  {
    std::cerr << "ERROR: cannot create "
              << outputPath << "\n";
    return;
  }

  out->cd();

  TNamed method(
    "analysis_method",
    "Same-event Lambda/Lambda-bar/K0S/Phi correlation histogram producer; "
    "Lambda spin observable follows arXiv:2506.05499; "
    "K0S/Phi combinations are control correlations; "
    "V0s use secondary-PCA momenta and Phi uses stored primary-PCA momenta.");
  method.Write();

  const TString identityTitle = TString::Format(
    "evt=%s, run=%s, daughter tracks=%s/%s",
    eventBranch, runBranch, trackId1Branch, trackId2Branch);
  TNamed identity(
    "identity_branches",
    identityTitle.Data());
  identity.Write();

  TNamed crossingDefinition(
    "crossing_definition",
    "A candidate has an internally resolved crossing only when "
    "crossing_valid1==1 AND crossing_valid2==1 AND crossing1==crossing2. "
    "Crossing value 0 is allowed when valid. "
    "sameCrossing/differentCrossing pair folders are filled only after both "
    "candidates satisfy this internal requirement.");
  crossingDefinition.Write();

  TNamed eventDefinition(
    "event_definition",
    "Candidates are grouped by input-file, run and evt. "
    "evt is the input frame/event sequence; multiple physical crossings inside "
    "the same evt are intentionally retained and separated by crossing QA.");
  eventDefinition.Write();

  TNamed rangeDef(
    "pair_range_definition",
    "short: |Delta y|<0.5 AND |Delta phi|<pi/3; "
    "long: (0.5<=|Delta y|<2) OR "
    "(pi/3<=|Delta phi|<pi)");
  rangeDef.Write();

  TParameter<int>(
      "mixed_event_enabled",
      doMixedEvent ? 1 : 0).Write();

  if (doMixedEvent)
  {
    const TString mixedDefTitle = TString::Format(
      "optional two-direction event mixing; relative "
      "|Delta pT|,|Delta phi|,|Delta y| matched within "
      "%.3f,%.3f,%.3f; each SE pair gives total ME weight 1",
      mixedMatchDeltaPt,
      mixedMatchDeltaPhi,
      mixedMatchDeltaY);
    TNamed mixedDef(
      "mixed_event_definition",
      mixedDefTitle.Data());
    mixedDef.Write();
  }

  TNamed controlAngles(
    "control_angle_definition",
    "For Lambda/anti-Lambda the analyzer is p/anti-p. "
    "For K0S controls it is pi+ (US) and for Phi controls K+ (US). "
    "The K0S/Phi mixed-species cos(theta*) histograms are control angular "
    "correlations, not a direct Lambda-Lambda spin observable.");
  controlAngles.Write();

  TParameter<double>(
      "alpha_Lambda",kAlphaLambda).Write();
  TParameter<double>(
      "alpha_AntiLambda",kAlphaAntiLambda).Write();

  HistogramStore hs;

  // ------------------------------------------------------------------------
  // Single-candidate QA
  // ------------------------------------------------------------------------
  TDirectory* candRoot =
      out->mkdir("candidates");

  std::vector<bool> signClasses{false};
  if (includeLikeSignTree)
    signClasses.push_back(true);

  for (bool like : signClasses)
  {
    TDirectory* signDir =
        candRoot->mkdir(like ? "LS" : "US");

    for (Species species : {
           Species::Lambda,
           Species::AntiLambda,
           Species::K0S,
           Species::Phi})
    {
      TDirectory* d =
          signDir->mkdir(speciesName(species));

      hs.candidate[
        HistogramStore::candidateKey(species,like)] =
          bookCandidateHists(d,species,like);
    }
  }

  // ------------------------------------------------------------------------
  // Pair histograms
  // ------------------------------------------------------------------------
  std::vector<PairSource> sources{
    PairSource::USUS
  };

  if (includeLikeSignTree)
  {
    sources.push_back(PairSource::USLS);
    sources.push_back(PairSource::LSLS);
  }

  // Do not even book the heavy mixed-event histogram hierarchy unless
  // explicitly requested.
  if (doMixedEvent)
    sources.push_back(PairSource::MixedEvent);

  const std::vector<PairType> pairTypes{
    PairType::LambdaLambda,
    PairType::LambdaAntiLambda,
    PairType::AntiLambdaAntiLambda,
    PairType::K0SK0S,
    PairType::K0SLambda,
    PairType::K0SAntiLambda,
    PairType::LambdaPhi,
    PairType::AntiLambdaPhi,
    PairType::K0SPhi
  };

  const std::vector<CrossingClass> crossingClasses{
    CrossingClass::All,
    CrossingClass::BothCandidatesDaughtersSame,
    CrossingClass::SameResolvedCrossing,
    CrossingClass::DifferentResolvedCrossing,
    CrossingClass::InternalMismatch,
    CrossingClass::InvalidCrossing
  };

  TDirectory* pairRoot =
      out->mkdir("pairs");

  for (PairSource source : sources)
  {
    TDirectory* sourceDir =
        pairRoot->mkdir(sourceName(source));

    for (PairType type : pairTypes)
    {
      TDirectory* typeDir =
          sourceDir->mkdir(pairTypeName(type));

      for (const std::string& range : {
             std::string("all"),
             std::string("short"),
             std::string("long")})
      {
        TDirectory* rangeDir =
            typeDir->mkdir(range.c_str());

        const std::vector<CrossingClass> classesToBook =
            source == PairSource::MixedEvent
              ? std::vector<CrossingClass>{CrossingClass::All}
              : crossingClasses;

        for (CrossingClass cc : classesToBook)
        {
          TDirectory* d =
              rangeDir->mkdir(
                  crossingClassName(cc));

          hs.pair[
            HistogramStore::pairKey(
              source,type,range,cc)] =
                bookPairHists(d,type);
        }
      }
    }
  }

  // ------------------------------------------------------------------------
  // Compact pair tree
  // ------------------------------------------------------------------------
  out->cd();

  PairSummaryRow summaryRow;
  TTree* summary = nullptr;

  if (writePairSummary)
  {
    summary =
      new TTree(
        "pairSummary",
        "Compact same-event pair observables and crossing QA");

    setupPairSummaryTree(
        summary,summaryRow);
  }

  TH1I* hEvents =
      new TH1I(
        "h_event_counts",
        ";category;events",
        4,0,4);

  hEvents->GetXaxis()->SetBinLabel(
      1,"input evt groups");
  hEvents->GetXaxis()->SetBinLabel(
      2,"evt with >=2 US candidates");
  hEvents->GetXaxis()->SetBinLabel(
      3,"evt with LS candidates");
  hEvents->GetXaxis()->SetBinLabel(
      4,"evt mixed");

  std::deque<EventRecord> mixPool;
  Long64_t processedEntries = 0;
  Long64_t totalEvents = 0;

  // ------------------------------------------------------------------------
  // Input processing. pairTree rows from one Fun4All event are consecutive,
  // and the event identity is stored in evt. We still key by file/run/evt so
  // a repeated evt number in a different input file cannot be merged.
  // ------------------------------------------------------------------------
  for (std::size_t fileIndex=0;
       fileIndex<files.size();
       ++fileIndex)
  {
    if (maxEntries >= 0 &&
        processedEntries >= maxEntries)
      break;

    std::unique_ptr<TFile> f(
        TFile::Open(
          files[fileIndex].c_str(),
          "READ"));

    if (!f || f->IsZombie())
    {
      std::cerr
        << "WARNING: cannot open "
        << files[fileIndex]
        << "\n";
      continue;
    }

    TTree* pairTree =
      dynamic_cast<TTree*>(
        f->Get(pairTreeName));

    if (!pairTree)
    {
      std::cerr
        << "WARNING: missing "
        << pairTreeName
        << " in "
        << files[fileIndex]
        << "\n";
      continue;
    }

    std::map<EventKey,EventRecord>
      events;

    addCandidatesFromTree(
        pairTree,
        false,
        static_cast<int>(fileIndex),
        events,
        hs,
        processedEntries,
        maxEntries,
        eventBranch,
        runBranch,
        trackId1Branch,
        trackId2Branch,
        requireFourDistinctTracks,
        applyCurrentSignedDeltaPhiCut,
        requireProtonHigherPt,
        beamX,
        beamY);

    if (includeLikeSignTree)
    {
      if (TTree* likeTree =
          dynamic_cast<TTree*>(
            f->Get(likeSignTreeName)))
      {
        addCandidatesFromTree(
            likeTree,
            true,
            static_cast<int>(fileIndex),
            events,
            hs,
            processedEntries,
            maxEntries,
            eventBranch,
            runBranch,
            trackId1Branch,
            trackId2Branch,
            requireFourDistinctTracks,
            false,
            false,
            beamX,
            beamY);
      }
      else
      {
        std::cerr
          << "WARNING: missing "
          << likeSignTreeName
          << " in "
          << files[fileIndex]
          << "\n";
      }
    }

    for (auto& kv : events)
    {
      EventRecord& event =
          kv.second;

      ++totalEvents;
      hEvents->Fill(0.5);

      const std::size_t nUS =
          event.usLambda.size() +
          event.usAntiLambda.size() +
          event.usK0S.size() +
          event.usPhi.size();

      const std::size_t nLS =
          event.lsLambda.size() +
          event.lsAntiLambda.size() +
          event.lsK0S.size() +
          event.lsPhi.size();

      if (nUS >= 2)
        hEvents->Fill(1.5);

      if (nLS > 0)
        hEvents->Fill(2.5);

      processSameEvent(
          event,
          hs,
          summary,
          summaryRow,
          requireFourDistinctTracks);

      if (doMixedEvent &&
          !mixPool.empty())
      {
        buildMixedEvent(
            event,
            mixPool,
            hs,
            requireFourDistinctTracks,
            mixedMatchDeltaPt,
            mixedMatchDeltaPhi,
            mixedMatchDeltaY);

        hEvents->Fill(3.5);
      }

      if (doMixedEvent)
      {
        mixPool.push_back(event);

        while (
          static_cast<int>(
            mixPool.size()) > mixDepth)
        {
          mixPool.pop_front();
        }
      }
    }

    std::cout
      << "Processed file "
      << (fileIndex+1)
      << "/"
      << files.size()
      << ": "
      << files[fileIndex]
      << ", cumulative tree entries="
      << processedEntries
      << ", evt groups="
      << totalEvents
      << "\n";
  }

  out->cd();

  TParameter<Long64_t>(
      "processed_tree_entries",
      processedEntries).Write();

  TParameter<Long64_t>(
      "processed_event_groups",
      totalEvents).Write();

  out->Write();
  out->Close();

  std::cout
    << "Wrote "
    << outputPath
    << "\n";

  std::cout
    << "Crossing split: all -> bothCandidatesDaughtersSame -> "
    << "sameCrossing/differentCrossing, with internalMismatch and "
    << "invalidCrossing kept as QA.\n";

  if (doMixedEvent)
  {
    std::cout
      << "Mixed-event histograms were enabled. "
      << "They can be normalized/divided downstream if desired.\n";
  }
  else
  {
    std::cout
      << "Mixed-event histograms were OFF (default). "
      << "No mixing CPU or mixed-event TH3 output was produced.\n";
  }
}
