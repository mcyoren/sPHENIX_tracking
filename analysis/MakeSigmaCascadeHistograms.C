// MakeSigmaCascadeHistograms.C
//
// One-pass QA/cut scan for TpcHadronReco:
//   sigma1385Tree : Sigma(1385) -> Lambda pi
//   cascadeTree   : Xi -> Lambda pi, Omega -> Lambda K
//
// NOTE: these trees already contain producer-level selected candidates.
// This macro can only apply tighter offline selections.

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TNamed.h>
#include <TParameter.h>
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
constexpr double kLambdaMass = 1.115683;
constexpr UInt_t kXi = 1U << 0;
constexpr UInt_t kOmega = 1U << 1;

bool has(TChain& t, const char* b) { return t.GetBranch(b) != nullptr; }

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

bool passSigma(const SigmaCut& c, double ml, double ldca, double ldira,
               double lr, double bpt, unsigned int bn, double bdca)
{
  if (c.dmLambda >= 0 && (!std::isfinite(ml) || std::abs(ml-kLambdaMass) > c.dmLambda)) return false;
  if (c.maxLambdaDca >= 0 && (!std::isfinite(ldca) || std::abs(ldca) > c.maxLambdaDca)) return false;
  if (c.minLambdaDira >= -1 && (!std::isfinite(ldira) || ldira < c.minLambdaDira)) return false;
  if (c.minLambdaR >= 0 && (!std::isfinite(lr) || lr < c.minLambdaR)) return false;
  if (c.minBachelorPt >= 0 && (!std::isfinite(bpt) || bpt < c.minBachelorPt)) return false;
  if (c.minBachelorN > 0 && static_cast<int>(bn) < c.minBachelorN) return false;
  if (c.maxBachelorDcaXY >= 0 && (!std::isfinite(bdca) || std::abs(bdca) > c.maxBachelorDcaXY)) return false;
  return true;
}

bool passCascade(const CascadeCut& c, double ml, double ldca, double bpt,
                 unsigned int bn, double cdca, double cr, double cdira,
                 double lflight, double ldira, double bdca)
{
  if (c.dmLambda >= 0 && (!std::isfinite(ml) || std::abs(ml-kLambdaMass) > c.dmLambda)) return false;
  if (c.maxLambdaDca >= 0 && (!std::isfinite(ldca) || std::abs(ldca) > c.maxLambdaDca)) return false;
  if (c.minBachelorPt >= 0 && (!std::isfinite(bpt) || bpt < c.minBachelorPt)) return false;
  if (c.minBachelorN > 0 && static_cast<int>(bn) < c.minBachelorN) return false;
  if (c.maxCascadeDca >= 0 && (!std::isfinite(cdca) || std::abs(cdca) > c.maxCascadeDca)) return false;
  if (c.minCascadeR >= 0 && (!std::isfinite(cr) || cr < c.minCascadeR)) return false;
  if (c.minCascadeDira >= -1 && (!std::isfinite(cdira) || cdira < c.minCascadeDira)) return false;
  if (c.minLambdaFlight >= 0 && (!std::isfinite(lflight) || lflight < c.minLambdaFlight)) return false;
  if (c.minLambdaDira >= -1 && (!std::isfinite(ldira) || ldira < c.minLambdaDira)) return false;
  if (c.minBachelorDcaXY >= 0 && (!std::isfinite(bdca) || std::abs(bdca) < c.minBachelorDcaXY)) return false;
  return true;
}

struct SigmaH
{
  TH1F *mass{}, *lamP{}, *lamM{}, *alamP{}, *alamM{};
  TH2F *pt{}, *lambdaMass{}, *lambdaR{}, *lambdaDca{}, *lambdaDira{}, *bachelorPt{}, *bachelorDca{};
};

SigmaH bookSigma(TDirectory* d, const std::string& tag)
{
  d->cd();
  SigmaH h;
  h.mass = new TH1F("h_mass_sigma1385",
      ("#Sigma(1385) #rightarrow #Lambda#pi ["+tag+"];m_{#Lambda#pi} [GeV/c^{2}];candidates").c_str(),
      300,1.25,1.55);
  h.lamP = new TH1F("h_mass_LambdaPiPlus","#Lambda#pi^{+};m [GeV/c^{2}];candidates",300,1.25,1.55);
  h.lamM = new TH1F("h_mass_LambdaPiMinus","#Lambda#pi^{-};m [GeV/c^{2}];candidates",300,1.25,1.55);
  h.alamP = new TH1F("h_mass_AntiLambdaPiPlus","#bar{#Lambda}#pi^{+};m [GeV/c^{2}];candidates",300,1.25,1.55);
  h.alamM = new TH1F("h_mass_AntiLambdaPiMinus","#bar{#Lambda}#pi^{-};m [GeV/c^{2}];candidates",300,1.25,1.55);
  h.pt = new TH2F("h_mass_vs_pt","mass vs p_{T};p_{T}^{#Lambda#pi} [GeV/c];m [GeV/c^{2}]",
      60,0,6,300,1.25,1.55);
  h.lambdaMass = new TH2F("h_mass_vs_lambda_mass","mass vs #Lambda mass;m_{#Lambda};m_{#Lambda#pi}",
      160,1.075,1.155,300,1.25,1.55);
  h.lambdaR = new TH2F("h_mass_vs_lambda_decayR","mass vs #Lambda decay radius;R_{#Lambda} [cm];m",
      100,0,50,300,1.25,1.55);
  h.lambdaDca = new TH2F("h_mass_vs_lambda_pairDCA","mass vs #Lambda pair DCA;DCA_{p#pi} [cm];m",
      100,0,3,300,1.25,1.55);
  h.lambdaDira = new TH2F("h_mass_vs_lambda_DIRA","mass vs #Lambda DIRA;DIRA_{#Lambda};m",
      120,-0.2,1.0,300,1.25,1.55);
  h.bachelorPt = new TH2F("h_mass_vs_bachelor_pt","mass vs bachelor p_{T};p_{T}^{#pi};m",
      60,0,3,300,1.25,1.55);
  h.bachelorDca = new TH2F("h_mass_vs_bachelor_dcaXY","mass vs bachelor |DCA_{xy}|;|DCA_{xy}| [cm];m",
      100,0,3,300,1.25,1.55);
  return h;
}

struct CascadeH
{
  TH1F *xi{}, *xim{}, *xip{}, *om{}, *omm{}, *omp{};
  TH2F *xiPt{}, *omPt{}, *xiDca{}, *omDca{}, *xiR{}, *omR{},
       *xiDira{}, *omDira{}, *xiLFlight{}, *omLFlight{},
       *xiLDira{}, *omLDira{}, *xiBDca{}, *omBDca{}, *xiLmass{}, *omLmass{};
};

CascadeH bookCascade(TDirectory* d)
{
  d->cd();
  CascadeH h;
  h.xi  = new TH1F("h_mass_xi","#Xi^{#mp};m_{#Lambda#pi} [GeV/c^{2}];candidates",300,1.25,1.40);
  h.xim = new TH1F("h_mass_xi_minus","#Xi^{-};m_{#Lambda#pi^{-}} [GeV/c^{2}];candidates",300,1.25,1.40);
  h.xip = new TH1F("h_mass_xi_plus","#bar{#Xi}^{+};m_{#bar{#Lambda}#pi^{+}} [GeV/c^{2}];candidates",300,1.25,1.40);
  h.om  = new TH1F("h_mass_omega","#Omega^{#mp};m_{#Lambda K} [GeV/c^{2}];candidates",300,1.60,1.75);
  h.omm = new TH1F("h_mass_omega_minus","#Omega^{-};m_{#Lambda K^{-}} [GeV/c^{2}];candidates",300,1.60,1.75);
  h.omp = new TH1F("h_mass_omega_plus","#bar{#Omega}^{+};m_{#bar{#Lambda}K^{+}} [GeV/c^{2}];candidates",300,1.60,1.75);

  h.xiPt = new TH2F("h_xi_mass_vs_pt",";p_{T}^{#Xi};m_{#Xi}",60,0,6,300,1.25,1.40);
  h.omPt = new TH2F("h_omega_mass_vs_pt",";p_{T}^{#Omega};m_{#Omega}",60,0,6,300,1.60,1.75);
  h.xiDca = new TH2F("h_xi_mass_vs_cascadeDCA",";DCA_{#Lambda,bach} [cm];m_{#Xi}",100,0,2,300,1.25,1.40);
  h.omDca = new TH2F("h_omega_mass_vs_cascadeDCA",";DCA_{#Lambda,bach} [cm];m_{#Omega}",100,0,2,300,1.60,1.75);
  h.xiR = new TH2F("h_xi_mass_vs_decayR",";R_{cascade} [cm];m_{#Xi}",100,0,50,300,1.25,1.40);
  h.omR = new TH2F("h_omega_mass_vs_decayR",";R_{cascade} [cm];m_{#Omega}",100,0,50,300,1.60,1.75);
  h.xiDira = new TH2F("h_xi_mass_vs_DIRA",";DIRA_{cascade};m_{#Xi}",100,0.70,1,300,1.25,1.40);
  h.omDira = new TH2F("h_omega_mass_vs_DIRA",";DIRA_{cascade};m_{#Omega}",100,0.70,1,300,1.60,1.75);
  h.xiLFlight = new TH2F("h_xi_mass_vs_lambdaFlight",";L_{#Lambda}^{cascade} [cm];m_{#Xi}",100,0,50,300,1.25,1.40);
  h.omLFlight = new TH2F("h_omega_mass_vs_lambdaFlight",";L_{#Lambda}^{cascade} [cm];m_{#Omega}",100,0,50,300,1.60,1.75);
  h.xiLDira = new TH2F("h_xi_mass_vs_lambdaDIRA",";DIRA_{#Lambda}^{cascade};m_{#Xi}",100,0.70,1,300,1.25,1.40);
  h.omLDira = new TH2F("h_omega_mass_vs_lambdaDIRA",";DIRA_{#Lambda}^{cascade};m_{#Omega}",100,0.70,1,300,1.60,1.75);
  h.xiBDca = new TH2F("h_xi_mass_vs_bachelorDCAxy",";|DCA_{xy}^{bach}| [cm];m_{#Xi}",100,0,5,300,1.25,1.40);
  h.omBDca = new TH2F("h_omega_mass_vs_bachelorDCAxy",";|DCA_{xy}^{bach}| [cm];m_{#Omega}",100,0,5,300,1.60,1.75);
  h.xiLmass = new TH2F("h_xi_mass_vs_lambdaMass",";m_{#Lambda};m_{#Xi}",160,1.075,1.155,300,1.25,1.40);
  h.omLmass = new TH2F("h_omega_mass_vs_lambdaMass",";m_{#Lambda};m_{#Omega}",160,1.075,1.155,300,1.60,1.75);
  return h;
}
}

void MakeSigmaCascadeHistograms(
    const char* inputDir=".",
    const char* filePattern="*.root",
    const char* outputDir="output",
    const char* outputName="sigma_cascade_histograms.root",
    const char* sigmaTreeName="sigma1385Tree",
    const char* cascadeTreeName="cascadeTree",
    const Long64_t maxEntries=-1,
    const double beamX=0.158,
    const double beamY=0.285)
{
  TH1::AddDirectory(kTRUE);
  const TString pattern = TString::Format("%s/%s",inputDir,filePattern);
  gSystem->mkdir(outputDir,kTRUE);
  const TString outpath = TString::Format("%s/%s",outputDir,outputName);

  std::unique_ptr<TFile> out(TFile::Open(outpath,"RECREATE"));
  if(!out || out->IsZombie()) { std::cerr<<"Cannot create "<<outpath<<std::endl; return; }

  TNamed("note","Trees are already producer-selected; this macro only tightens cuts.").Write();
  TParameter<double>("beam_x_cm",beamX).Write();
  TParameter<double>("beam_y_cm",beamY).Write();

  // ---------------- Sigma(1385) ----------------
  TChain s(sigmaTreeName);
  const int nsfiles=s.Add(pattern);
  if(nsfiles>0 && s.GetEntries()>0 && has(s,"mass_sigma1385"))
  {
    const std::vector<std::string> req={
      "is_antilambda","sigma_charge","bachelor_ntpc_clusters","lambda_mass",
      "lambda_pair_dca","lambda_dira","lambda_decay_x","lambda_decay_y",
      "bachelor_pt","bachelor_dca_xy","sigma_pt","mass_sigma1385"};
    for(const auto& b:req) if(!has(s,b.c_str())) { std::cerr<<"Missing Sigma branch "<<b<<std::endl; return; }

    Int_t anti=0,charge=0; UInt_t ncl=0;
    Float_t lm=0,ldca=0,ldira=0,lx=0,ly=0,bpt=0,bdca=0,spt=0,m=0;
    s.SetBranchAddress("is_antilambda",&anti);
    s.SetBranchAddress("sigma_charge",&charge);
    s.SetBranchAddress("bachelor_ntpc_clusters",&ncl);
    s.SetBranchAddress("lambda_mass",&lm);
    s.SetBranchAddress("lambda_pair_dca",&ldca);
    s.SetBranchAddress("lambda_dira",&ldira);
    s.SetBranchAddress("lambda_decay_x",&lx);
    s.SetBranchAddress("lambda_decay_y",&ly);
    s.SetBranchAddress("bachelor_pt",&bpt);
    s.SetBranchAddress("bachelor_dca_xy",&bdca);
    s.SetBranchAddress("sigma_pt",&spt);
    s.SetBranchAddress("mass_sigma1385",&m);

    const std::vector<SigmaCut> cuts={
      {"cut00_raw","producer output",-1,-1,-2,-1,-1,0,-1},
      {"cut01_loose","dmL<40MeV, LpairDCA<3, RL>1, bach pT>0.2, nTPC>=20",.040,3,-2,1,.20,20,3},
      {"cut02","dmL<30MeV, LpairDCA<2, LDIRA>0.75, RL>1",.030,2,.75,1,.20,20,3},
      {"cut03_baseline","dmL<25MeV, LpairDCA<1.5, LDIRA>0.85, RL>1.5, nTPC>=25",.025,1.5,.85,1.5,.20,25,2.5},
      {"cut04","dmL<20MeV, LpairDCA<1, LDIRA>0.90, RL>2, nTPC>=30",.020,1,.90,2,.20,30,2},
      {"cut05","dmL<15MeV, LpairDCA<0.8, LDIRA>0.95, RL>2, bach pT>0.25",.015,.8,.95,2,.25,30,1.5},
      {"cut06","dmL<12MeV, LpairDCA<0.6, LDIRA>0.97, RL>2.5",.012,.6,.97,2.5,.25,30,1},
      {"cut07","dmL<10MeV, LpairDCA<0.5, LDIRA>0.98, RL>3, bach pT>0.30",.010,.5,.98,3,.30,35,1},
      {"cut08","dmL<8MeV, LpairDCA<0.4, LDIRA>0.99, RL>3",.008,.4,.99,3,.30,35,.7},
      {"cut09_tight","dmL<6MeV, LpairDCA<0.3, LDIRA>0.995, RL>3",.006,.3,.995,3,.30,35,.5}
    };

    TDirectory* top=out->mkdir("sigma1385");
    top->cd();
    TH1I* flow=new TH1I("h_cutflow","Sigma cumulative cut flow;cut;candidates",cuts.size(),0,cuts.size());
    std::map<std::string,SigmaH> hs;
    for(size_t i=0;i<cuts.size();++i){
      flow->GetXaxis()->SetBinLabel(i+1,cuts[i].name.c_str());
      auto* d=top->mkdir(cuts[i].name.c_str()); d->cd();
      TNamed("selection",cuts[i].desc.c_str()).Write();
      hs[cuts[i].name]=bookSigma(d,cuts[i].name);
    }

    const Long64_t n=maxEntries<0?s.GetEntries():std::min(maxEntries,s.GetEntries());
    for(Long64_t ie=0;ie<n;++ie){
      s.GetEntry(ie);
      const double lr=std::hypot(lx-beamX,ly-beamY);
      for(size_t ic=0;ic<cuts.size();++ic){
        const auto& c=cuts[ic];
        if(!passSigma(c,lm,ldca,ldira,lr,bpt,ncl,bdca)) continue;
        flow->AddBinContent(ic+1);
        auto& h=hs[c.name];
        h.mass->Fill(m); h.pt->Fill(spt,m); h.lambdaMass->Fill(lm,m); h.lambdaR->Fill(lr,m);
        h.lambdaDca->Fill(std::abs(ldca),m); h.lambdaDira->Fill(ldira,m);
        h.bachelorPt->Fill(bpt,m); h.bachelorDca->Fill(std::abs(bdca),m);
        if(!anti && charge>0) h.lamP->Fill(m);
        else if(!anti && charge<0) h.lamM->Fill(m);
        else if(anti && charge>0) h.alamP->Fill(m);
        else if(anti && charge<0) h.alamM->Fill(m);
      }
    }
    std::cout<<"Sigma files "<<nsfiles<<", entries "<<s.GetEntries()<<std::endl;
  } else std::cout<<"No usable "<<sigmaTreeName<<" in "<<pattern<<std::endl;

  // ---------------- Xi / Omega ----------------
  TChain c(cascadeTreeName);
  const int ncfiles=c.Add(pattern);
  if(ncfiles>0 && c.GetEntries()>0 && has(c,"candidate_mask"))
  {
    const std::vector<std::string> req={
      "is_antilambda","cascade_charge","candidate_mask","bachelor_ntpc_clusters",
      "lambda_mass","lambda_pair_dca","cascade_pair_dca","cascade_decay_radius",
      "cascade_dira","lambda_flight_from_cascade","lambda_dira_from_cascade",
      "bachelor_pt","bachelor_dca_xy","cascade_pt","mass_xi","mass_omega"};
    for(const auto& b:req) if(!has(c,b.c_str())) { std::cerr<<"Missing cascade branch "<<b<<std::endl; return; }

    Int_t anti=0,charge=0; UInt_t mask=0,ncl=0;
    Float_t lm=0,ldca=0,cdca=0,cr=0,cdira=0,lf=0,ldira=0,bpt=0,bdca=0,cpt=0,mxi=0,mom=0;
    c.SetBranchAddress("is_antilambda",&anti);
    c.SetBranchAddress("cascade_charge",&charge);
    c.SetBranchAddress("candidate_mask",&mask);
    c.SetBranchAddress("bachelor_ntpc_clusters",&ncl);
    c.SetBranchAddress("lambda_mass",&lm);
    c.SetBranchAddress("lambda_pair_dca",&ldca);
    c.SetBranchAddress("cascade_pair_dca",&cdca);
    c.SetBranchAddress("cascade_decay_radius",&cr);
    c.SetBranchAddress("cascade_dira",&cdira);
    c.SetBranchAddress("lambda_flight_from_cascade",&lf);
    c.SetBranchAddress("lambda_dira_from_cascade",&ldira);
    c.SetBranchAddress("bachelor_pt",&bpt);
    c.SetBranchAddress("bachelor_dca_xy",&bdca);
    c.SetBranchAddress("cascade_pt",&cpt);
    c.SetBranchAddress("mass_xi",&mxi);
    c.SetBranchAddress("mass_omega",&mom);

    const std::vector<CascadeCut> cuts={
      {"cut00_raw","producer output",-1,-1,-1,0,-1,-1,-2,-1,-2,-1},
      {"cut01_producer_like","dmL<40MeV, LDCA<3, bach pT>0.2, nTPC>=20",.040,3,.20,20,2,.5,.75,.2,.75,-1},
      {"cut02_loose","dmL<30MeV, cascadeDCA<1.5, R>0.7, DIRA>0.80",.030,2.5,.20,20,1.5,.7,.80,.3,.80,-1},
      {"cut03_baseline","dmL<25MeV, cascadeDCA<1, R>1, DIRA>0.85, Llambda>0.5",.025,2,.20,25,1,1,.85,.5,.85,.02},
      {"cut04","dmL<20MeV, cascadeDCA<0.8, DIRA>0.90, nTPC>=30",.020,1.5,.20,30,.8,1,.90,.7,.90,.03},
      {"cut05","LDCA<1, cascadeDCA<0.6, R>1.5, DIRA>0.93, bach pT>0.25",.020,1,.25,30,.6,1.5,.93,1,.93,.05},
      {"cut06","dmL<15MeV, cascadeDCA<0.5, DIRA>0.95",.015,.8,.25,30,.5,1.5,.95,1,.95,.05},
      {"cut07","cascadeDCA<0.4, R>2, DIRA>0.97, Llambda>1.5, nTPC>=35",.015,.6,.30,35,.4,2,.97,1.5,.97,.08},
      {"cut08","dmL<12MeV, cascadeDCA<0.3, DIRA>0.98, Llambda>2",.012,.5,.30,35,.3,2,.98,2,.98,.10},
      {"cut09_tight","dmL<10MeV, cascadeDCA<0.2, R>2.5, DIRA>0.99",.010,.4,.30,35,.2,2.5,.99,2,.99,.10}
    };

    TDirectory* top=out->mkdir("cascade"); top->cd();
    TH1I* fxi=new TH1I("h_cutflow_xi","Xi cumulative cut flow;cut;candidates",cuts.size(),0,cuts.size());
    TH1I* fom=new TH1I("h_cutflow_omega","Omega cumulative cut flow;cut;candidates",cuts.size(),0,cuts.size());
    std::map<std::string,CascadeH> hc;
    for(size_t i=0;i<cuts.size();++i){
      fxi->GetXaxis()->SetBinLabel(i+1,cuts[i].name.c_str());
      fom->GetXaxis()->SetBinLabel(i+1,cuts[i].name.c_str());
      auto* d=top->mkdir(cuts[i].name.c_str()); d->cd();
      TNamed("selection",cuts[i].desc.c_str()).Write();
      hc[cuts[i].name]=bookCascade(d);
    }

    const Long64_t n=maxEntries<0?c.GetEntries():std::min(maxEntries,c.GetEntries());
    for(Long64_t ie=0;ie<n;++ie){
      c.GetEntry(ie);
      const bool xi=mask&kXi, om=mask&kOmega;
      for(size_t ic=0;ic<cuts.size();++ic){
        const auto& cut=cuts[ic];
        if(!passCascade(cut,lm,ldca,bpt,ncl,cdca,cr,cdira,lf,ldira,bdca)) continue;
        auto& h=hc[cut.name];
        if(xi){
          fxi->AddBinContent(ic+1);
          h.xi->Fill(mxi); h.xiPt->Fill(cpt,mxi); h.xiDca->Fill(std::abs(cdca),mxi);
          h.xiR->Fill(cr,mxi); h.xiDira->Fill(cdira,mxi); h.xiLFlight->Fill(lf,mxi);
          h.xiLDira->Fill(ldira,mxi); h.xiBDca->Fill(std::abs(bdca),mxi); h.xiLmass->Fill(lm,mxi);
          if(!anti && charge<0) h.xim->Fill(mxi); else if(anti && charge>0) h.xip->Fill(mxi);
        }
        if(om){
          fom->AddBinContent(ic+1);
          h.om->Fill(mom); h.omPt->Fill(cpt,mom); h.omDca->Fill(std::abs(cdca),mom);
          h.omR->Fill(cr,mom); h.omDira->Fill(cdira,mom); h.omLFlight->Fill(lf,mom);
          h.omLDira->Fill(ldira,mom); h.omBDca->Fill(std::abs(bdca),mom); h.omLmass->Fill(lm,mom);
          if(!anti && charge<0) h.omm->Fill(mom); else if(anti && charge>0) h.omp->Fill(mom);
        }
      }
    }
    std::cout<<"Cascade files "<<ncfiles<<", entries "<<c.GetEntries()<<std::endl;
  } else std::cout<<"No usable "<<cascadeTreeName<<" in "<<pattern<<std::endl;

  out->Write();
  out->Close();
  std::cout<<"Wrote "<<outpath<<std::endl;
}
