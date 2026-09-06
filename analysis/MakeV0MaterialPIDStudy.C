// MakeV0MaterialPIDStudy.C
//
// Focused TPC-only test of whether the radial V0 background structures are
// detector-material hadronic interactions, especially
//
//      pi + A -> pi + p + X
//
// using the existing TpcHadronReco pairTree.
//
// The macro compares:
//
//   1) pi+ pi- PID pairs (K0S-like control)
//   2) p+ pi- PID pairs (secondary-proton / material hypothesis)
//   3) pi+ pbar- PID pairs (charge-conjugate control)
//
// while also scanning topology variables:
//   opening angle, |Delta phi|, |Delta eta|, pT asymmetry,
//   pair DCA, daughter DCA product, DIRA, Armenteros alpha/qT.
//
// dE/dx PID:
//   - ALEPH/ALICE center curve from the supplied sPHENIX study.
//   - Nsigma widths use the universal width parameterization
//       sigma_U(x) = sqrt(a^2 + b^2/x^c)
//     with a=0.07690, b=0.16914, c=0.40000.
//   - K, p, d use the per-species multiplicative width corrections from
//     the presentation, frozen outside their measured momentum ranges.
//   - The presentation's final pion width is a kernel smoother requiring
//     the individual pion width loci, which are not tabulated in the PDF.
//     Therefore this macro uses sigma_U itself for the pion width. This is
//     explicitly marked as the pion-width fallback.
//
// Strict PID follows the presentation logic:
//   accept |Nsigma_X| < pidNSigma
//   and veto all other implemented species with |Nsigma_Y| > pidNSigma.
//
// Material coordinates are DETECTOR-centered:
//   R = hypot(pca_x,pca_y), phi = atan2(pca_y,pca_x)
//
// Example:
// root -l -b -q \
// 'MakeV0MaterialPIDStudy.C("/path/to/v0","v0_pp_*.root","v0_material_pid.root")'
//

#include <TBranch.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TMath.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  constexpr double kPi = TMath::Pi();

  constexpr double kMassPi = 0.13957039;
  constexpr double kMassK  = 0.493677;
  constexpr double kMassP  = 0.938272088;
  constexpr double kMassD  = 1.87561294257;

  // ALEPH center-curve parameters.
  constexpr double kP1 = 0.1819;
  constexpr double kP2 = 10.8175;
  constexpr double kP3 = 1.1674e-5;
  constexpr double kP4 = 3.0199;
  constexpr double kP5 = 5.0;

  // Universal width in z = ln(dE/dx_meas)-ln(dE/dx_exp).
  constexpr double kWidthA = 0.07690;
  constexpr double kWidthB = 0.16914;
  constexpr double kWidthC = 0.40000;

  enum class Species
  {
    Pion = 0,
    Kaon = 1,
    Proton = 2,
    Deuteron = 3
  };

  double massOf(const Species s)
  {
    switch (s)
    {
      case Species::Pion: return kMassPi;
      case Species::Kaon: return kMassK;
      case Species::Proton: return kMassP;
      case Species::Deuteron: return kMassD;
    }
    return kMassPi;
  }

  double p3(const double px,const double py,const double pz)
  {
    return std::sqrt(px*px+py*py+pz*pz);
  }

  double pt(const double px,const double py)
  {
    return std::hypot(px,py);
  }

  double wrapPhi(const double x)
  {
    return std::atan2(std::sin(x),std::cos(x));
  }

  double eta(const double px,const double py,const double pz)
  {
    const double pT=std::hypot(px,py);
    if(!(pT>0.0)) return std::numeric_limits<double>::quiet_NaN();
    return std::asinh(pz/pT);
  }

  double openingAngle(const double px1,const double py1,const double pz1,
                      const double px2,const double py2,const double pz2)
  {
    const double p1=p3(px1,py1,pz1);
    const double p2=p3(px2,py2,pz2);
    if(!(p1>0.0) || !(p2>0.0))
      return std::numeric_limits<double>::quiet_NaN();

    double c=(px1*px2+py1*py2+pz1*pz2)/(p1*p2);
    c=std::clamp(c,-1.0,1.0);
    return std::acos(c);
  }

  double asymmetry(const double a,const double b)
  {
    return (a+b)>0.0
      ? std::abs(a-b)/(a+b)
      : std::numeric_limits<double>::quiet_NaN();
  }

  double invariantMass(const double px1,const double py1,const double pz1,
                       const double m1,
                       const double px2,const double py2,const double pz2,
                       const double m2)
  {
    const double pp1=px1*px1+py1*py1+pz1*pz1;
    const double pp2=px2*px2+py2*py2+pz2*pz2;

    const double e1=std::sqrt(pp1+m1*m1);
    const double e2=std::sqrt(pp2+m2*m2);

    const double px=px1+px2;
    const double py=py1+py2;
    const double pz=pz1+pz2;

    const double m2pair=(e1+e2)*(e1+e2)-(px*px+py*py+pz*pz);
    return m2pair>0.0 ? std::sqrt(m2pair) : 0.0;
  }

  // --------------------------------------------------------------------------
  // dE/dx center and width parameterizations
  // --------------------------------------------------------------------------

  double alephExpectedDedx(const double betaGamma)
  {
    if(!(betaGamma>0.0))
      return std::numeric_limits<double>::quiet_NaN();

    const double beta=betaGamma/std::sqrt(1.0+betaGamma*betaGamma);
    const double betaP4=std::pow(beta,kP4);
    const double arg=kP3+std::pow(betaGamma,-kP5);

    if(!(betaP4>0.0) || !(arg>0.0))
      return std::numeric_limits<double>::quiet_NaN();

    const double x=(kP1/betaP4)*(kP2-betaP4-std::log(arg));
    return x>0.0 ? x : std::numeric_limits<double>::quiet_NaN();
  }

  double expectedDedx(const double p,const Species species)
  {
    if(!(p>0.0))
      return std::numeric_limits<double>::quiet_NaN();

    return alephExpectedDedx(p/massOf(species));
  }

  double expectedLnDedx(const double p,const Species species)
  {
    const double x=expectedDedx(p,species);
    return x>0.0 ? std::log(x)
                 : std::numeric_limits<double>::quiet_NaN();
  }

  double sigmaUniversal(const double expectedDedxValue)
  {
    if(!(expectedDedxValue>0.0))
      return std::numeric_limits<double>::quiet_NaN();

    return std::sqrt(
      kWidthA*kWidthA +
      kWidthB*kWidthB/std::pow(expectedDedxValue,kWidthC));
  }

  double clampedLogMomentum(const double p,
                            const double pLo,
                            const double pHi)
  {
    const double pc=std::clamp(p,pLo,pHi);
    return std::log(pc);
  }

  double speciesWidthRatio(const double p,const Species species)
  {
    // The PDF gives no tabulated pion kernel-smoothed width points, so use
    // the universal shape itself for pion: r_pi = 1.
    if(species==Species::Pion) return 1.0;

    if(species==Species::Kaon)
    {
      constexpr double pLo=0.23;
      constexpr double pHi=0.53;
      const double t=clampedLogMomentum(p,pLo,pHi);
      return 0.910 + 0.628*(t-std::log(pLo));
    }

    if(species==Species::Proton)
    {
      constexpr double pLo=0.20;
      constexpr double pHi=0.93;
      const double t=clampedLogMomentum(p,pLo,pHi);
      return 1.077 - 0.037*(t-std::log(pLo));
    }

    if(species==Species::Deuteron)
    {
      constexpr double pLo=0.29;
      constexpr double pHi=1.73;
      const double t=clampedLogMomentum(p,pLo,pHi);
      return 0.834 + 0.164*(t-std::log(pLo));
    }

    return 1.0;
  }

  double sigmaSpecies(const double p,const Species species)
  {
    const double x=expectedDedx(p,species);
    const double su=sigmaUniversal(x);
    if(!std::isfinite(su)) return su;
    return speciesWidthRatio(p,species)*su;
  }

  struct PID
  {
    bool valid{false};
    double lnDedx{std::numeric_limits<double>::quiet_NaN()};

    std::array<double,4> nSigma{{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()
    }};

    bool strictPion{false};
    bool strictProton{false};
    bool loosePion{false};
    bool looseProton{false};
  };

  std::size_t indexOf(const Species species)
  {
    return static_cast<std::size_t>(species);
  }

  bool strictSpeciesPID(const PID& pid,
                        const Species target,
                        const double nSigmaCut)
  {
    if(!pid.valid) return false;

    for(const Species s:
        {Species::Pion,Species::Kaon,Species::Proton,Species::Deuteron})
    {
      const double ns=pid.nSigma[indexOf(s)];
      if(!std::isfinite(ns)) return false;

      if(s==target)
      {
        if(std::abs(ns)>=nSigmaCut) return false;
      }
      else
      {
        if(std::abs(ns)<=nSigmaCut) return false;
      }
    }
    return true;
  }

  PID makePID(const double p,
              const double dedx,
              const double dedxScale,
              const double nSigmaCut)
  {
    PID out;
    if(!(p>0.0) || !(dedx>0.0) || !(dedxScale>0.0))
      return out;

    out.lnDedx=std::log(dedx*dedxScale);

    for(const Species s:
        {Species::Pion,Species::Kaon,Species::Proton,Species::Deuteron})
    {
      const double mu=expectedLnDedx(p,s);
      const double sigma=sigmaSpecies(p,s);
      if(!std::isfinite(mu) || !(sigma>0.0))
        return out;

      out.nSigma[indexOf(s)]=(out.lnDedx-mu)/sigma;
    }

    out.valid=true;

    out.loosePion=
      std::abs(out.nSigma[indexOf(Species::Pion)])<nSigmaCut;
    out.looseProton=
      std::abs(out.nSigma[indexOf(Species::Proton)])<nSigmaCut;

    out.strictPion=strictSpeciesPID(out,Species::Pion,nSigmaCut);
    out.strictProton=strictSpeciesPID(out,Species::Proton,nSigmaCut);

    return out;
  }

  struct Pair
  {
    // ordered by electric charge
    double pxPos{},pyPos{},pzPos{};
    double pxNeg{},pyNeg{},pzNeg{};
    double dedxPos{},dedxNeg{};

    double pPos{},pNeg{},ptPos{},ptNeg{};

    PID pidPos,pidNeg;

    double massPiPi{};
    double massPPiMinus{};   // p+ + pi-
    double massPiPlusPbar{}; // pi+ + pbar-

    bool piPiPID{false};
    bool pPiMinusPID{false};
    bool piPbarPID{false};

    double x{},y{},z{},r{},phi{};
    double opening{},absDphi{},absDeta{},absDtheta{},phiEtaRatio{};
    double ptAsym{},pAsym{};
    double pairDca{},dira{},lproj{},alpha{},qT{};
    double dcaPos{},dcaNeg{},minAbsDca{},signedDcaProduct{},absDcaProduct{};
    double maxQuality{};
    int minNpoints{};
  };

  struct Topology
  {
    std::string name;
    std::string description;

    double maxPairDca{-1.0};
    double minAbsDca{-1.0};
    double minAbsDcaProduct{-1.0};
    double minDira{-2.0};

    double minOpening{-1.0};
    double minAbsDphi{-1.0};
    double maxAbsDeta{-1.0};
    double minPhiEtaRatio{-1.0};

    double minPtAsym{-1.0};
    double maxPtAsym{-1.0};

    double minAbsAlpha{-1.0};
    double maxAbsAlpha{-1.0};
    double maxQT{-1.0};

    int minNpoints{0};
    double maxQuality{-1.0};
  };

  bool passTopology(const Pair& p,const Topology& c)
  {
    if(c.maxPairDca>=0.0 && p.pairDca>c.maxPairDca) return false;
    if(c.minAbsDca>=0.0 && p.minAbsDca<c.minAbsDca) return false;
    if(c.minAbsDcaProduct>=0.0 && p.absDcaProduct<c.minAbsDcaProduct) return false;
    if(c.minDira>=-1.0 && p.dira<c.minDira) return false;

    if(c.minOpening>=0.0 && p.opening<c.minOpening) return false;
    if(c.minAbsDphi>=0.0 && p.absDphi<c.minAbsDphi) return false;
    if(c.maxAbsDeta>=0.0 && p.absDeta>c.maxAbsDeta) return false;
    if(c.minPhiEtaRatio>=0.0 && p.phiEtaRatio<c.minPhiEtaRatio) return false;

    if(c.minPtAsym>=0.0 && p.ptAsym<c.minPtAsym) return false;
    if(c.maxPtAsym>=0.0 && p.ptAsym>c.maxPtAsym) return false;

    if(c.minAbsAlpha>=0.0 && std::abs(p.alpha)<c.minAbsAlpha) return false;
    if(c.maxAbsAlpha>=0.0 && std::abs(p.alpha)>c.maxAbsAlpha) return false;
    if(c.maxQT>=0.0 && p.qT>c.maxQT) return false;

    if(c.minNpoints>0 && p.minNpoints<c.minNpoints) return false;
    if(c.maxQuality>=0.0 && p.maxQuality>c.maxQuality) return false;

    return true;
  }

  struct Maps
  {
    TH2D* xy{};
    TH2D* rphi{};
    TH2D* rz{};
    TH2D* phiz{};
    TH3F* rphiz{};
  };

  Maps bookMaps(TDirectory* d,
                const std::string& prefix,
                const std::string& title)
  {
    d->cd();
    Maps h;

    h.xy=new TH2D(
      (prefix+"_xy").c_str(),
      (title+";x_{sec} [cm];y_{sec} [cm]").c_str(),
      240,-30,30,240,-30,30);

    h.rphi=new TH2D(
      (prefix+"_rphi").c_str(),
      (title+";#phi_{sec};R_{sec} [cm]").c_str(),
      180,-kPi,kPi,140,0,35);

    h.rz=new TH2D(
      (prefix+"_rz").c_str(),
      (title+";z_{sec} [cm];R_{sec} [cm]").c_str(),
      240,-120,120,140,0,35);

    h.phiz=new TH2D(
      (prefix+"_phiz").c_str(),
      (title+";z_{sec} [cm];#phi_{sec}").c_str(),
      240,-120,120,180,-kPi,kPi);

    h.rphiz=new TH3F(
      (prefix+"_rphiz").c_str(),
      (title+";R_{sec} [cm];#phi_{sec};z_{sec} [cm]").c_str(),
      70,0,35,72,-kPi,kPi,120,-120,120);

    return h;
  }

  void fillMaps(Maps& h,const Pair& p)
  {
    h.xy->Fill(p.x,p.y);
    h.rphi->Fill(p.phi,p.r);
    h.rz->Fill(p.z,p.r);
    h.phiz->Fill(p.z,p.phi);
    h.rphiz->Fill(p.r,p.phi,p.z);
  }

  struct CategoryHist
  {
    TH1D* massPiPi{};
    TH1D* massPPiMinus{};
    TH1D* massPiPlusPbar{};

    TH2D* piPiVsPPi{};
    TH2D* massPiPiVsR{};
    TH2D* massPiPiVsOpening{};
    TH2D* massPiPiVsPtAsym{};
    TH2D* massPiPiVsDphi{};
    TH2D* massPiPiVsDeta{};
    TH2D* dphiVsDeta{};
    TH2D* massPiPiVsPairDca{};
    TH2D* massPiPiVsDira{};
    TH2D* massPiPiVsDcaProduct{};
    TH2D* massPiPiVsAlpha{};
    TH2D* massPiPiVsQT{};

    Maps allMaps;
    Maps k0WindowMaps;
  };

  CategoryHist bookCategory(TDirectory* d,
                            const std::string& title)
  {
    d->cd();
    CategoryHist h;

    h.massPiPi=new TH1D(
      "h_mass_pipi",
      (title+";m_{#pi#pi} [GeV/c^{2}];pairs").c_str(),
      400,0.40,0.60);

    h.massPPiMinus=new TH1D(
      "h_mass_pPiMinus",
      (title+";m_{p^{+}#pi^{-}} [GeV/c^{2}];pairs").c_str(),
      500,1.05,2.05);

    h.massPiPlusPbar=new TH1D(
      "h_mass_PiPlusPbar",
      (title+";m_{#pi^{+}#bar{p}} [GeV/c^{2}];pairs").c_str(),
      500,1.05,2.05);

    h.piPiVsPPi=new TH2D(
      "h_mass_pipi_vs_pPi",
      (title+";m_{#pi#pi} [GeV/c^{2}];m_{p#pi} [GeV/c^{2}]").c_str(),
      240,0.42,0.58,300,1.05,1.80);

    h.massPiPiVsR=new TH2D(
      "h_mass_pipi_vs_R",
      (title+";R_{sec} [cm];m_{#pi#pi} [GeV/c^{2}]").c_str(),
      140,0,35,320,0.42,0.58);

    h.massPiPiVsOpening=new TH2D(
      "h_mass_pipi_vs_opening",
      (title+";opening angle [rad];m_{#pi#pi}").c_str(),
      180,0,kPi,320,0.42,0.58);

    h.massPiPiVsPtAsym=new TH2D(
      "h_mass_pipi_vs_ptAsym",
      (title+";A_{pT};m_{#pi#pi}").c_str(),
      100,0,1,320,0.42,0.58);

    h.massPiPiVsDphi=new TH2D(
      "h_mass_pipi_vs_absDphi",
      (title+";|#Delta#phi| [rad];m_{#pi#pi}").c_str(),
      180,0,kPi,320,0.42,0.58);

    h.massPiPiVsDeta=new TH2D(
      "h_mass_pipi_vs_absDeta",
      (title+";|#Delta#eta|;m_{#pi#pi}").c_str(),
      160,0,4,320,0.42,0.58);

    h.dphiVsDeta=new TH2D(
      "h_absDphi_vs_absDeta",
      (title+";|#Delta#eta|;|#Delta#phi| [rad]").c_str(),
      160,0,4,180,0,kPi);

    h.massPiPiVsPairDca=new TH2D(
      "h_mass_pipi_vs_pairDCA",
      (title+";pair DCA [cm];m_{#pi#pi}").c_str(),
      160,0,4,320,0.42,0.58);

    h.massPiPiVsDira=new TH2D(
      "h_mass_pipi_vs_DIRA",
      (title+";DIRA;m_{#pi#pi}").c_str(),
      200,-1,1,320,0.42,0.58);

    h.massPiPiVsDcaProduct=new TH2D(
      "h_mass_pipi_vs_absDCAProduct",
      (title+";|DCA_{xy,+}DCA_{xy,-}| [cm^{2}];m_{#pi#pi}").c_str(),
      160,0,4,320,0.42,0.58);

    h.massPiPiVsAlpha=new TH2D(
      "h_mass_pipi_vs_alpha",
      (title+";#alpha;m_{#pi#pi}").c_str(),
      160,-1.2,1.2,320,0.42,0.58);

    h.massPiPiVsQT=new TH2D(
      "h_mass_pipi_vs_qT",
      (title+";q_{T} [GeV/c];m_{#pi#pi}").c_str(),
      160,0,0.4,320,0.42,0.58);

    h.allMaps=bookMaps(d->mkdir("all_mass_maps"),"h_all",title+" all masses");
    h.k0WindowMaps=bookMaps(
      d->mkdir("k0_window_maps"),"h_k0window",
      title+" with 0.488 < m_{#pi#pi} < 0.508");

    return h;
  }

  void fillCategory(CategoryHist& h,
                    const Pair& p,
                    const bool inK0Window)
  {
    h.massPiPi->Fill(p.massPiPi);
    h.massPPiMinus->Fill(p.massPPiMinus);
    h.massPiPlusPbar->Fill(p.massPiPlusPbar);

    // For the physical p+pi- hypothesis use m(p+ pi-).
    h.piPiVsPPi->Fill(p.massPiPi,p.massPPiMinus);

    h.massPiPiVsR->Fill(p.r,p.massPiPi);
    h.massPiPiVsOpening->Fill(p.opening,p.massPiPi);
    h.massPiPiVsPtAsym->Fill(p.ptAsym,p.massPiPi);
    h.massPiPiVsDphi->Fill(p.absDphi,p.massPiPi);
    if(std::isfinite(p.absDeta))
    {
      h.massPiPiVsDeta->Fill(p.absDeta,p.massPiPi);
      h.dphiVsDeta->Fill(p.absDeta,p.absDphi);
    }
    h.massPiPiVsPairDca->Fill(p.pairDca,p.massPiPi);
    h.massPiPiVsDira->Fill(p.dira,p.massPiPi);
    h.massPiPiVsDcaProduct->Fill(p.absDcaProduct,p.massPiPi);
    h.massPiPiVsAlpha->Fill(p.alpha,p.massPiPi);
    h.massPiPiVsQT->Fill(p.qT,p.massPiPi);

    fillMaps(h.allMaps,p);
    if(inK0Window) fillMaps(h.k0WindowMaps,p);
  }

  bool hasBranch(TChain& chain,const char* name)
  {
    return chain.GetBranch(name)!=nullptr;
  }
}

void MakeV0MaterialPIDStudy(
  const char* inputDir=".",
  const char* filePattern="*.root",
  const char* outputFile="v0_material_pid.root",
  const char* treeName="pairTree",
  Long64_t maxEntries=-1,

  // Stored dE/dx -> equalized dE/dx global scale.
  double dedxScale=0.015,

  // Nsigma threshold used for strict PID.
  double pidNSigma=2.0)
{
  const TString pattern=TString::Format("%s/%s",inputDir,filePattern);

  TChain chain(treeName);
  const int nFiles=chain.Add(pattern);
  if(nFiles<=0)
  {
    std::cerr<<"ERROR: no files matched "<<pattern<<std::endl;
    return;
  }

  const std::vector<std::string> required={
    "pca_x","pca_y","pca_z",
    "px1","py1","pz1","px2","py2","pz2",
    "charge1","charge2","dedx_1","dedx_2",
    "pairDCA","cosThetaReco","Lproj","alpha","qT",
    "dca_xy1","dca_xy2",
    "quality1","quality2","npoints1","npoints2"
  };

  for(const auto& b:required)
  {
    if(!hasBranch(chain,b.c_str()))
    {
      std::cerr<<"ERROR: missing branch "<<b<<std::endl;
      return;
    }
  }

  Float_t pca_x=0,pca_y=0,pca_z=0;
  Float_t px1=0,py1=0,pz1=0,px2=0,py2=0,pz2=0;
  Float_t charge1=0,charge2=0,dedx1=0,dedx2=0;
  Float_t pairDCA=0,cosThetaReco=0,Lproj=0,alpha=0,qT=0;
  Float_t dca_xy1=0,dca_xy2=0,quality1=0,quality2=0;
  Short_t npoints1=0,npoints2=0;

  chain.SetBranchAddress("pca_x",&pca_x);
  chain.SetBranchAddress("pca_y",&pca_y);
  chain.SetBranchAddress("pca_z",&pca_z);

  chain.SetBranchAddress("px1",&px1);
  chain.SetBranchAddress("py1",&py1);
  chain.SetBranchAddress("pz1",&pz1);
  chain.SetBranchAddress("px2",&px2);
  chain.SetBranchAddress("py2",&py2);
  chain.SetBranchAddress("pz2",&pz2);

  chain.SetBranchAddress("charge1",&charge1);
  chain.SetBranchAddress("charge2",&charge2);
  chain.SetBranchAddress("dedx_1",&dedx1);
  chain.SetBranchAddress("dedx_2",&dedx2);

  chain.SetBranchAddress("pairDCA",&pairDCA);
  chain.SetBranchAddress("cosThetaReco",&cosThetaReco);
  chain.SetBranchAddress("Lproj",&Lproj);
  chain.SetBranchAddress("alpha",&alpha);
  chain.SetBranchAddress("qT",&qT);

  chain.SetBranchAddress("dca_xy1",&dca_xy1);
  chain.SetBranchAddress("dca_xy2",&dca_xy2);
  chain.SetBranchAddress("quality1",&quality1);
  chain.SetBranchAddress("quality2",&quality2);
  chain.SetBranchAddress("npoints1",&npoints1);
  chain.SetBranchAddress("npoints2",&npoints2);

  // Topology scans.  No mass or PID is built into these.
  const std::vector<Topology> topologies={
    {
      "topo00_all",
      "all unlike-sign pairs"
    },

    {
      "topo01_good_vertex",
      "pairDCA<1, min|DCAxy|>0.03, npoints>=20",
      1.0,0.03,-1.0,-2.0,
      -1.0,-1.0,-1.0,-1.0,
      -1.0,-1.0,
      -1.0,-1.0,-1.0,
      20,-1.0
    },

    {
      "topo02_phi_open",
      "good vertex + opening>0.08, |dphi|>0.04, |dphi|/(|deta|+1e-3)>1",
      1.0,0.03,-1.0,-2.0,
      0.08,0.04,-1.0,1.0,
      -1.0,-1.0,
      -1.0,-1.0,-1.0,
      20,-1.0
    },

    {
      "topo03_material_loose",
      "pairDCA<0.7, min|DCAxy|>0.05, opening>0.10, |dphi|>0.05, phi/eta>1, pT asym>0.30",
      0.7,0.05,-1.0,-2.0,
      0.10,0.05,-1.0,1.0,
      0.30,-1.0,
      -1.0,-1.0,-1.0,
      20,-1.0
    },

    {
      "topo04_material_dca_product",
      "material loose + |DCAxy+*DCAxy-|>0.01 cm2",
      0.7,0.05,0.01,-2.0,
      0.10,0.05,-1.0,1.0,
      0.30,-1.0,
      -1.0,-1.0,-1.0,
      20,-1.0
    },

    {
      "topo05_k0s",
      "K0S-like: pairDCA<0.5, min|DCAxy|>0.05, DIRA>0.95, |alpha|<0.8, qT<0.22, npoints>=25",
      0.5,0.05,-1.0,0.95,
      -1.0,-1.0,-1.0,-1.0,
      -1.0,0.85,
      -1.0,0.8,0.22,
      25,-1.0
    },

    {
      "topo06_high_asym",
      "hadronic-interaction stress test: pairDCA<1, min|DCAxy|>0.05, pT asym>0.55, |alpha|>0.45",
      1.0,0.05,-1.0,-2.0,
      -1.0,-1.0,-1.0,-1.0,
      0.55,-1.0,
      0.45,-1.0,-1.0,
      20,-1.0
    }
  };

  std::unique_ptr<TFile> out(TFile::Open(outputFile,"RECREATE"));
  if(!out || out->IsZombie())
  {
    std::cerr<<"ERROR: cannot create "<<outputFile<<std::endl;
    return;
  }

  out->cd();

  TNamed("input",pattern.Data()).Write();
  TNamed(
    "pid_note",
    "Strict PID uses |Nsigma_X|<cut and |Nsigma_Y|>cut for pi,K,p,d. "
    "Pion width uses universal sigma_U because the PDF does not tabulate the "
    "individual pion width loci required to reproduce its final kernel smoother.").Write();

  TParameter<double>("dedxScale",dedxScale).Write();
  TParameter<double>("pidNSigma",pidNSigma).Write();

  TParameter<double>("ALEPH_P1",kP1).Write();
  TParameter<double>("ALEPH_P2",kP2).Write();
  TParameter<double>("ALEPH_P3",kP3).Write();
  TParameter<double>("ALEPH_P4",kP4).Write();
  TParameter<double>("ALEPH_P5",kP5).Write();

  TParameter<double>("width_a",kWidthA).Write();
  TParameter<double>("width_b",kWidthB).Write();
  TParameter<double>("width_c",kWidthC).Write();

  // Global PID QA: this is the first thing to inspect to validate dedxScale.
  TH2D* hLnDedxVsP=new TH2D(
    "h_lnDedx_vs_p",
    "all daughters;p [GeV/c];ln(dE/dx)",
    220,0,5.5,220,0,5.5);

  TH2D* hNsPiVsP=new TH2D(
    "h_NsigmaPi_vs_p",
    "all daughters;p [GeV/c];N#sigma_{#pi}",
    220,0,5.5,240,-12,12);

  TH2D* hNsKVsP=new TH2D(
    "h_NsigmaK_vs_p",
    "all daughters;p [GeV/c];N#sigma_{K}",
    220,0,5.5,240,-12,12);

  TH2D* hNsPVsP=new TH2D(
    "h_NsigmaP_vs_p",
    "all daughters;p [GeV/c];N#sigma_{p}",
    220,0,5.5,240,-12,12);

  TH2D* hNsDVsP=new TH2D(
    "h_NsigmaD_vs_p",
    "all daughters;p [GeV/c];N#sigma_{d}",
    220,0,5.5,240,-12,12);

  TH2D* hNsPiVsPPos=new TH2D(
    "h_NsigmaPi_vs_NsigmaP_positive",
    "positive daughters;N#sigma_{#pi};N#sigma_{p}",
    240,-12,12,240,-12,12);

  TH2D* hNsPiVsPNeg=new TH2D(
    "h_NsigmaPi_vs_NsigmaP_negative",
    "negative daughters;N#sigma_{#pi};N#sigma_{p}",
    240,-12,12,240,-12,12);

  // One directory per topology; inside each are the three PID hypotheses.
  struct TopologyOutput
  {
    CategoryHist all;
    CategoryHist piPi;
    CategoryHist pPiMinus;
    CategoryHist piPbar;
  };

  std::map<std::string,TopologyOutput> hTopo;

  for(const auto& topo:topologies)
  {
    TDirectory* dTopo=out->mkdir(topo.name.c_str());
    dTopo->cd();
    TNamed("topology",topo.description.c_str()).Write();

    TopologyOutput h;
    h.all=bookCategory(dTopo->mkdir("pid00_all"),"all unlike-sign");
    h.piPi=bookCategory(
      dTopo->mkdir("pid01_piPlus_piMinus"),
      "strict #pi^{+}#pi^{-} PID");

    h.pPiMinus=bookCategory(
      dTopo->mkdir("pid02_pPlus_piMinus"),
      "strict p^{+}#pi^{-} PID");

    h.piPbar=bookCategory(
      dTopo->mkdir("pid03_piPlus_pbarMinus"),
      "strict #pi^{+}#bar{p}^{-} PID control");

    hTopo[topo.name]=h;
  }

  TH1I* hCounts=new TH1I(
    "h_category_counts",
    "PID category counts;category;pairs",
    4,0,4);

  hCounts->GetXaxis()->SetBinLabel(1,"all unlike");
  hCounts->GetXaxis()->SetBinLabel(2,"pi+ pi-");
  hCounts->GetXaxis()->SetBinLabel(3,"p+ pi-");
  hCounts->GetXaxis()->SetBinLabel(4,"pi+ pbar-");

  // Compact derived tree for interactive follow-up.
  TTree derived("materialPIDTree","derived V0/material PID variables");

  Float_t o_mPiPi=0,o_mPPi=0,o_mPiPbar=0;
  Float_t o_x=0,o_y=0,o_z=0,o_r=0,o_phi=0;
  Float_t o_open=0,o_dphi=0,o_deta=0,o_ptAsym=0,o_pAsym=0;
  Float_t o_pairDca=0,o_dira=0,o_dcaProd=0,o_alpha=0,o_qT=0;
  Float_t o_pPos=0,o_pNeg=0,o_dedxPos=0,o_dedxNeg=0;
  Float_t o_nsPiPos=0,o_nsPPos=0,o_nsPiNeg=0,o_nsPNeg=0;
  Int_t o_piPiPID=0,o_pPiPID=0,o_piPbarPID=0;
  ULong64_t o_topologyMask=0;

  derived.Branch("massPiPi",&o_mPiPi,"massPiPi/F");
  derived.Branch("massPPiMinus",&o_mPPi,"massPPiMinus/F");
  derived.Branch("massPiPlusPbar",&o_mPiPbar,"massPiPlusPbar/F");

  derived.Branch("x",&o_x,"x/F");
  derived.Branch("y",&o_y,"y/F");
  derived.Branch("z",&o_z,"z/F");
  derived.Branch("r",&o_r,"r/F");
  derived.Branch("phi",&o_phi,"phi/F");

  derived.Branch("opening",&o_open,"opening/F");
  derived.Branch("absDphi",&o_dphi,"absDphi/F");
  derived.Branch("absDeta",&o_deta,"absDeta/F");
  derived.Branch("ptAsym",&o_ptAsym,"ptAsym/F");
  derived.Branch("pAsym",&o_pAsym,"pAsym/F");

  derived.Branch("pairDCA",&o_pairDca,"pairDCA/F");
  derived.Branch("DIRA",&o_dira,"DIRA/F");
  derived.Branch("absDCAProduct",&o_dcaProd,"absDCAProduct/F");
  derived.Branch("alpha",&o_alpha,"alpha/F");
  derived.Branch("qT",&o_qT,"qT/F");

  derived.Branch("pPositive",&o_pPos,"pPositive/F");
  derived.Branch("pNegative",&o_pNeg,"pNegative/F");
  derived.Branch("dedxPositive",&o_dedxPos,"dedxPositive/F");
  derived.Branch("dedxNegative",&o_dedxNeg,"dedxNegative/F");

  derived.Branch("NsigmaPiPositive",&o_nsPiPos,"NsigmaPiPositive/F");
  derived.Branch("NsigmaPPositive",&o_nsPPos,"NsigmaPPositive/F");
  derived.Branch("NsigmaPiNegative",&o_nsPiNeg,"NsigmaPiNegative/F");
  derived.Branch("NsigmaPNegative",&o_nsPNeg,"NsigmaPNegative/F");

  derived.Branch("isPiPiPID",&o_piPiPID,"isPiPiPID/I");
  derived.Branch("isPPiMinusPID",&o_pPiPID,"isPPiMinusPID/I");
  derived.Branch("isPiPlusPbarPID",&o_piPbarPID,"isPiPlusPbarPID/I");
  derived.Branch("topologyMask",&o_topologyMask,"topologyMask/l");

  const Long64_t nEntries=
    maxEntries<0 ? chain.GetEntries()
                 : std::min(maxEntries,chain.GetEntries());

  Long64_t nUnlike=0;

  for(Long64_t ie=0;ie<nEntries;++ie)
  {
    chain.GetEntry(ie);

    // This study is charge-specific and starts from unlike-sign pairs.
    if(charge1*charge2>=0) continue;
    ++nUnlike;

    Pair p;

    const bool oneIsPositive=charge1>0;

    if(oneIsPositive)
    {
      p.pxPos=px1; p.pyPos=py1; p.pzPos=pz1; p.dedxPos=dedx1;
      p.pxNeg=px2; p.pyNeg=py2; p.pzNeg=pz2; p.dedxNeg=dedx2;
      p.dcaPos=dca_xy1; p.dcaNeg=dca_xy2;
    }
    else
    {
      p.pxPos=px2; p.pyPos=py2; p.pzPos=pz2; p.dedxPos=dedx2;
      p.pxNeg=px1; p.pyNeg=py1; p.pzNeg=pz1; p.dedxNeg=dedx1;
      p.dcaPos=dca_xy2; p.dcaNeg=dca_xy1;
    }

    p.pPos=p3(p.pxPos,p.pyPos,p.pzPos);
    p.pNeg=p3(p.pxNeg,p.pyNeg,p.pzNeg);
    p.ptPos=pt(p.pxPos,p.pyPos);
    p.ptNeg=pt(p.pxNeg,p.pyNeg);

    p.pidPos=makePID(p.pPos,p.dedxPos,dedxScale,pidNSigma);
    p.pidNeg=makePID(p.pNeg,p.dedxNeg,dedxScale,pidNSigma);

    // Charge-specific hypotheses.
    p.piPiPID=p.pidPos.strictPion && p.pidNeg.strictPion;
    p.pPiMinusPID=p.pidPos.strictProton && p.pidNeg.strictPion;
    p.piPbarPID=p.pidPos.strictPion && p.pidNeg.strictProton;

    p.massPiPi=invariantMass(
      p.pxPos,p.pyPos,p.pzPos,kMassPi,
      p.pxNeg,p.pyNeg,p.pzNeg,kMassPi);

    p.massPPiMinus=invariantMass(
      p.pxPos,p.pyPos,p.pzPos,kMassP,
      p.pxNeg,p.pyNeg,p.pzNeg,kMassPi);

    p.massPiPlusPbar=invariantMass(
      p.pxPos,p.pyPos,p.pzPos,kMassPi,
      p.pxNeg,p.pyNeg,p.pzNeg,kMassP);

    p.x=pca_x;
    p.y=pca_y;
    p.z=pca_z;
    p.r=std::hypot(p.x,p.y);
    p.phi=std::atan2(p.y,p.x);

    p.opening=openingAngle(
      p.pxPos,p.pyPos,p.pzPos,
      p.pxNeg,p.pyNeg,p.pzNeg);

    p.absDphi=std::abs(
      wrapPhi(
        std::atan2(p.pyPos,p.pxPos)-
        std::atan2(p.pyNeg,p.pxNeg)));

    const double etaPos=eta(p.pxPos,p.pyPos,p.pzPos);
    const double etaNeg=eta(p.pxNeg,p.pyNeg,p.pzNeg);

    p.absDeta=
      std::isfinite(etaPos) && std::isfinite(etaNeg)
        ? std::abs(etaPos-etaNeg)
        : std::numeric_limits<double>::quiet_NaN();

    const double thetaPos=std::atan2(p.ptPos,p.pzPos);
    const double thetaNeg=std::atan2(p.ptNeg,p.pzNeg);
    p.absDtheta=std::abs(thetaPos-thetaNeg);

    p.phiEtaRatio=
      std::isfinite(p.absDeta)
        ? p.absDphi/(p.absDeta+1e-3)
        : 0.0;

    p.ptAsym=asymmetry(p.ptPos,p.ptNeg);
    p.pAsym=asymmetry(p.pPos,p.pNeg);

    p.pairDca=std::abs(pairDCA);
    p.dira=cosThetaReco;
    p.lproj=Lproj;
    p.alpha=alpha;
    p.qT=qT;

    p.dcaPos=std::abs(p.dcaPos);
    p.dcaNeg=std::abs(p.dcaNeg);
    p.minAbsDca=std::min(p.dcaPos,p.dcaNeg);

    // Preserve sign separately through input before abs would be better, but
    // this study primarily uses the magnitude of the product.
    p.absDcaProduct=p.dcaPos*p.dcaNeg;

    p.maxQuality=std::max(double(quality1),double(quality2));
    p.minNpoints=std::min(int(npoints1),int(npoints2));

    // Global dE/dx QA.
    for(const auto& item:
        {std::pair<double,PID>{p.pPos,p.pidPos},
         std::pair<double,PID>{p.pNeg,p.pidNeg}})
    {
      const double mom=item.first;
      const PID& pid=item.second;
      if(!pid.valid) continue;

      hLnDedxVsP->Fill(mom,pid.lnDedx);
      hNsPiVsP->Fill(mom,pid.nSigma[indexOf(Species::Pion)]);
      hNsKVsP->Fill(mom,pid.nSigma[indexOf(Species::Kaon)]);
      hNsPVsP->Fill(mom,pid.nSigma[indexOf(Species::Proton)]);
      hNsDVsP->Fill(mom,pid.nSigma[indexOf(Species::Deuteron)]);
    }

    if(p.pidPos.valid)
      hNsPiVsPPos->Fill(
        p.pidPos.nSigma[indexOf(Species::Pion)],
        p.pidPos.nSigma[indexOf(Species::Proton)]);

    if(p.pidNeg.valid)
      hNsPiVsPNeg->Fill(
        p.pidNeg.nSigma[indexOf(Species::Pion)],
        p.pidNeg.nSigma[indexOf(Species::Proton)]);

    hCounts->AddBinContent(1);
    if(p.piPiPID) hCounts->AddBinContent(2);
    if(p.pPiMinusPID) hCounts->AddBinContent(3);
    if(p.piPbarPID) hCounts->AddBinContent(4);

    const bool inK0Window=(p.massPiPi>=0.488 && p.massPiPi<0.508);

    o_topologyMask=0;

    for(std::size_t it=0;it<topologies.size();++it)
    {
      const auto& topo=topologies[it];
      if(!passTopology(p,topo)) continue;

      if(it<64) o_topologyMask|=(ULong64_t(1)<<it);

      auto& h=hTopo[topo.name];

      fillCategory(h.all,p,inK0Window);
      if(p.piPiPID) fillCategory(h.piPi,p,inK0Window);
      if(p.pPiMinusPID) fillCategory(h.pPiMinus,p,inK0Window);
      if(p.piPbarPID) fillCategory(h.piPbar,p,inK0Window);
    }

    o_mPiPi=p.massPiPi;
    o_mPPi=p.massPPiMinus;
    o_mPiPbar=p.massPiPlusPbar;

    o_x=p.x; o_y=p.y; o_z=p.z; o_r=p.r; o_phi=p.phi;
    o_open=p.opening; o_dphi=p.absDphi; o_deta=p.absDeta;
    o_ptAsym=p.ptAsym; o_pAsym=p.pAsym;
    o_pairDca=p.pairDca; o_dira=p.dira;
    o_dcaProd=p.absDcaProduct; o_alpha=p.alpha; o_qT=p.qT;

    o_pPos=p.pPos; o_pNeg=p.pNeg;
    o_dedxPos=p.dedxPos; o_dedxNeg=p.dedxNeg;

    o_nsPiPos=p.pidPos.valid
      ? p.pidPos.nSigma[indexOf(Species::Pion)] : -999.f;
    o_nsPPos=p.pidPos.valid
      ? p.pidPos.nSigma[indexOf(Species::Proton)] : -999.f;
    o_nsPiNeg=p.pidNeg.valid
      ? p.pidNeg.nSigma[indexOf(Species::Pion)] : -999.f;
    o_nsPNeg=p.pidNeg.valid
      ? p.pidNeg.nSigma[indexOf(Species::Proton)] : -999.f;

    o_piPiPID=p.piPiPID ? 1:0;
    o_pPiPID=p.pPiMinusPID ? 1:0;
    o_piPbarPID=p.piPbarPID ? 1:0;

    derived.Fill();
  }

  out->cd();
  //derived.Write();
  hCounts->Write();

  TNamed(
    "interpretation_hint",
    "If detector-material knockout dominates, compare p+pi- against pi+pbar-: "
    "the p+pi- category should show stronger fixed-R structures, especially "
    "after material-like topology cuts. True K0S should be concentrated in "
    "the pi+pi- PID category and remain symmetric in charge.").Write();

  out->Write();
  out->Close();

  std::cout
    <<"MakeV0MaterialPIDStudy finished\n"
    <<"  files: "<<nFiles<<"\n"
    <<"  entries read: "<<nEntries<<"\n"
    <<"  unlike pairs: "<<nUnlike<<"\n"
    <<"  output: "<<outputFile<<std::endl;
}
