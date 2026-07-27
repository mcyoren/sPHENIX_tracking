// TestK0sAlternatingConstraint.C
// Low-pT K0S -> pi+pi- alternating mass/pointing and pair-DCA constraint test.

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace
{
  constexpr double kPionMass = 0.13957039;
  constexpr double kKshortMass = 0.497611;

  struct Vec3
  {
    double x = 0.;
    double y = 0.;
    double z = 0.;
  };

  Vec3 add(const Vec3& a, const Vec3& b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
  Vec3 sub(const Vec3& a, const Vec3& b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
  Vec3 mul(const Vec3& a, double s) { return {s*a.x,s*a.y,s*a.z}; }
  double dot(const Vec3& a, const Vec3& b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
  double mag2(const Vec3& a) { return dot(a,a); }
  double mag(const Vec3& a) { return std::sqrt(mag2(a)); }
  double pt(const Vec3& a) { return std::hypot(a.x,a.y); }
  bool finite(const Vec3& a) { return std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(a.z); }

  Vec3 unit(const Vec3& a)
  {
    const double m = mag(a);
    if (!(m>0.))
    {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      return {nan,nan,nan};
    }
    return mul(a,1./m);
  }

  Vec3 rotateZ(const Vec3& p, double dphi)
  {
    const double c=std::cos(dphi), s=std::sin(dphi);
    return {c*p.x-s*p.y,s*p.x+c*p.y,p.z};
  }

  double massPipi(const Vec3& p1,const Vec3& p2)
  {
    const double e1=std::sqrt(mag2(p1)+kPionMass*kPionMass);
    const double e2=std::sqrt(mag2(p2)+kPionMass*kPionMass);
    const double m2=(e1+e2)*(e1+e2)-mag2(add(p1,p2));
    return m2>0. ? std::sqrt(m2) : 0.;
  }

  double pointing(const Vec3& p,const Vec3& flight)
  {
    const double den=mag(p)*mag(flight);
    return den>0. ? std::clamp(dot(p,flight)/den,-1.,1.) : -2.;
  }

  double openingAngle(const Vec3& p1,const Vec3& p2)
  {
    const double den=mag(p1)*mag(p2);
    if (!(den>0.)) return std::numeric_limits<double>::quiet_NaN();
    return std::acos(std::clamp(dot(p1,p2)/den,-1.,1.));
  }

  double lineLineDca(const Vec3& x1,const Vec3& p1,const Vec3& x2,const Vec3& p2)
  {
    const Vec3 u1=unit(p1), u2=unit(p2);
    if (!finite(u1)||!finite(u2)) return std::numeric_limits<double>::quiet_NaN();
    const Vec3 w=sub(x1,x2);
    const double a=dot(u1,u1), b=dot(u1,u2), c=dot(u2,u2);
    const double d=dot(u1,w), e=dot(u2,w);
    const double den=a*c-b*b;
    double s=0.,t=0.;
    if (std::abs(den)>1e-12)
    {
      s=(b*e-c*d)/den;
      t=(a*e-b*d)/den;
    }
    else if (c>0.) t=e/c;
    return mag(sub(add(x1,mul(u1,s)),add(x2,mul(u2,t))));
  }

  struct ScaleSolution
  {
    bool valid=false;
    double s1=1.,s2=1.;
  };

  ScaleSolution solveScales(const Vec3& p1,const Vec3& p2,const Vec3& flight,
                            double minScale,double maxScale)
  {
    ScaleSolution result;
    const double c1=p1.x*flight.y-p1.y*flight.x;
    const double c2=p2.x*flight.y-p2.y*flight.x;
    if (std::abs(c2)<1e-12) return result;
    const double ratio=-c1/c2;
    if (!(ratio>0.)||!std::isfinite(ratio)) return result;

    double best=std::numeric_limits<double>::max();
    for (int i=0;i<=1000;++i)
    {
      const double s1=minScale+(maxScale-minScale)*double(i)/1000.;
      const double s2=ratio*s1;
      if (s2<minScale||s2>maxScale) continue;
      const double dm=massPipi(mul(p1,s1),mul(p2,s2))-kKshortMass;
      const double penalty=1e-4*((s1-1.)*(s1-1.)+(s2-1.)*(s2-1.));
      const double obj=dm*dm+penalty;
      if (obj<best)
      {
        best=obj;
        result.s1=s1;
        result.s2=s2;
        result.valid=true;
      }
    }
    return result;
  }

  struct RotationSolution
  {
    bool valid=false;
    double dphi1=0.,dphi2=0.,dca=0.;
  };

  RotationSolution solveRotations(const Vec3& x1,const Vec3& x2,
                                  const Vec3& p1,const Vec3& p2,
                                  double maxAbsRotation)
  {
    RotationSolution result;
    double c1=0.,c2=0.,half=maxAbsRotation;
    for (int level=0;level<5;++level)
    {
      double best=std::numeric_limits<double>::max();
      double b1=c1,b2=c2,bdca=0.;
      for (int i=0;i<=40;++i)
      {
        const double d1=std::clamp(c1-half+2.*half*double(i)/40.,-maxAbsRotation,maxAbsRotation);
        for (int j=0;j<=40;++j)
        {
          const double d2=std::clamp(c2-half+2.*half*double(j)/40.,-maxAbsRotation,maxAbsRotation);
          const double dca=lineLineDca(x1,rotateZ(p1,d1),x2,rotateZ(p2,d2));
          if (!std::isfinite(dca)) continue;
          const double reg=0.05*(d1*d1+d2*d2)/(maxAbsRotation*maxAbsRotation);
          const double obj=(dca/0.05)*(dca/0.05)+reg;
          if (obj<best) { best=obj;b1=d1;b2=d2;bdca=dca; }
        }
      }
      if (!std::isfinite(best)) return result;
      c1=b1;c2=b2;half*=0.2;
      result={true,c1,c2,bdca};
    }
    return result;
  }
}

void TestK0sAlternatingConstraint(
    const char* inputDir=".",
    const char* filePattern="*.root",
    const char* outputName="k0s_alternating_constraint.root",
    const char* treeName="pairTree",
    const double beamX=0.158,
    const double beamY=0.285,
    const double beamZ=0.0,
    const double daughterPtMin=0.20,
    const double daughterPtMax=1.50,
    const double kshortPtMax=2.0,
    const double openingAngleMin=0.50,
    const double massMin=0.47,
    const double massMax=0.53,
    const double diraMin=0.90,
    const double pairDcaMax=1.0,
    const double maxAbsRotation=0.05,
    const double minScale=0.60,
    const double maxScale=1.40,
    const int iterations=3,
    const Long64_t maxEntries=-1)
{
  const TString pattern=TString::Format("%s/%s",inputDir,filePattern);
  TChain chain(treeName);
  if (chain.Add(pattern)<=0)
  {
    std::cerr<<"ERROR: no files matched "<<pattern<<std::endl;
    return;
  }

  const std::vector<const char*> required={
    "candidate_mask","charge1","charge2","mass_Kshort",
    "pca_x","pca_y","pca_z","pca1_x","pca1_y","pca1_z",
    "pca2_x","pca2_y","pca2_z","px1","py1","pz1","px2","py2","pz2",
    "pairDCA","quality1","quality2","npoints1","npoints2"};
  for (const char* name:required)
  {
    if (!chain.GetBranch(name))
    {
      std::cerr<<"ERROR: missing branch "<<name<<std::endl;
      return;
    }
  }

  UInt_t candidateMask=0;
  Float_t charge1=0,charge2=0,massK=0;
  Float_t pcaX=0,pcaY=0,pcaZ=0,pca1X=0,pca1Y=0,pca1Z=0,pca2X=0,pca2Y=0,pca2Z=0;
  Float_t px1=0,py1=0,pz1=0,px2=0,py2=0,pz2=0,pairDca=0,quality1=0,quality2=0;
  Short_t npoints1=0,npoints2=0;

#define B(name,var) chain.SetBranchAddress(name,&var)
  B("candidate_mask",candidateMask); B("charge1",charge1); B("charge2",charge2);
  B("mass_Kshort",massK); B("pca_x",pcaX); B("pca_y",pcaY); B("pca_z",pcaZ);
  B("pca1_x",pca1X); B("pca1_y",pca1Y); B("pca1_z",pca1Z);
  B("pca2_x",pca2X); B("pca2_y",pca2Y); B("pca2_z",pca2Z);
  B("px1",px1); B("py1",py1); B("pz1",pz1); B("px2",px2); B("py2",py2); B("pz2",pz2);
  B("pairDCA",pairDca); B("quality1",quality1); B("quality2",quality2);
  B("npoints1",npoints1); B("npoints2",npoints2);
#undef B

  const TString outputDir=gSystem->DirName(outputName);
  if (outputDir.Length()>0) gSystem->mkdir(outputDir,true);
  std::unique_ptr<TFile> output(TFile::Open(outputName,"RECREATE"));
  if (!output||output->IsZombie()) return;

  TH1D hMassBefore("h_mass_before","mass before;m_{#pi#pi} [GeV/c^{2}];candidates",240,0.44,0.56);
  TH1D hMassAfter("h_mass_after","mass after;m_{#pi#pi} [GeV/c^{2}];candidates",240,0.44,0.56);
  TH1D hDiraBefore("h_dira_before","DIRA before;DIRA;candidates",200,0.8,1.0);
  TH1D hDiraAfter("h_dira_after","DIRA after;DIRA;candidates",200,0.8,1.0);
  TH1D hDcaBefore("h_pair_dca_before","pair DCA before;pair DCA [cm];candidates",200,0,2);
  TH1D hDcaAfter("h_pair_dca_after","pair DCA after;pair DCA [cm];candidates",200,0,2);
  TH1D hS1("h_scale1","daughter 1 scale;s_{1};candidates",200,0.6,1.4);
  TH1D hS2("h_scale2","daughter 2 scale;s_{2};candidates",200,0.6,1.4);
  TH1D hR1("h_delta_phi1","daughter 1 rotation;#delta#phi_{1};candidates",200,-0.05,0.05);
  TH1D hR2("h_delta_phi2","daughter 2 rotation;#delta#phi_{2};candidates",200,-0.05,0.05);
  TH2D hScales("h_scale2_vs_scale1","scale correlation;s_{1};s_{2}",160,0.6,1.4,160,0.6,1.4);
  TH2D hRotations("h_delta_phi2_vs_delta_phi1","rotation correlation;#delta#phi_{1};#delta#phi_{2}",160,-0.05,0.05,160,-0.05,0.05);

  TTree result("constraintTree","alternating KShort constraint results");
  Float_t outMassBefore=0,outMassAfter=0,outDiraBefore=0,outDiraAfter=0;
  Float_t outDcaBefore=0,outDcaAfter=0,outScale1=1,outScale2=1,outRot1=0,outRot2=0;
  Float_t outPt1=0,outPt2=0,outOpening=0;
  Int_t outValid=0,outIterations=0;
#define O(name,var,type) result.Branch(name,&var,name "/" type)
  O("mass_before",outMassBefore,"F"); O("mass_after",outMassAfter,"F");
  O("dira_before",outDiraBefore,"F"); O("dira_after",outDiraAfter,"F");
  O("pair_dca_before",outDcaBefore,"F"); O("pair_dca_after",outDcaAfter,"F");
  O("scale1",outScale1,"F"); O("scale2",outScale2,"F");
  O("delta_phi1",outRot1,"F"); O("delta_phi2",outRot2,"F");
  O("pt1",outPt1,"F"); O("pt2",outPt2,"F"); O("opening_angle",outOpening,"F");
  O("valid",outValid,"I"); O("iterations",outIterations,"I");
#undef O

  const Long64_t nTotal=chain.GetEntries();
  const Long64_t nRun=maxEntries>0?std::min(maxEntries,nTotal):nTotal;
  Long64_t selected=0,solved=0;

  for (Long64_t entry=0;entry<nRun;++entry)
  {
    chain.GetEntry(entry);
    if (entry%100000==0) std::cout<<"Processing "<<entry<<" / "<<nRun<<std::endl;
    if ((candidateMask&1U)==0U||charge1*charge2>=0) continue;

    const Vec3 p1{px1,py1,pz1},p2{px2,py2,pz2};
    const Vec3 x1{pca1X,pca1Y,pca1Z},x2{pca2X,pca2Y,pca2Z};
    const Vec3 flight{sub(Vec3{pcaX,pcaY,pcaZ},Vec3{beamX,beamY,beamZ})};
    const double dpt1=pt(p1),dpt2=pt(p2),kpt=pt(add(p1,p2));
    const double angle=openingAngle(p1,p2),dir=pointing(add(p1,p2),flight);

    if (!finite(p1)||!finite(p2)||!std::isfinite(angle)||
        dpt1<daughterPtMin||dpt2<daughterPtMin||dpt1>daughterPtMax||dpt2>daughterPtMax||
        kpt>kshortPtMax||angle<openingAngleMin||massK<massMin||massK>massMax||
        dir<diraMin||std::abs(pairDca)>pairDcaMax||quality1>15||quality2>15||
        npoints1<=30||npoints2<=30) continue;

    ++selected;
    Vec3 c1=p1,c2=p2;
    double totalS1=1.,totalS2=1.,totalR1=0.,totalR2=0.;
    bool valid=true;
    int done=0;

    for (int iter=0;iter<iterations;++iter)
    {
      const ScaleSolution scales=solveScales(c1,c2,flight,minScale,maxScale);
      if (!scales.valid) { valid=false; break; }
      c1=mul(c1,scales.s1); c2=mul(c2,scales.s2);
      totalS1*=scales.s1; totalS2*=scales.s2;

      const RotationSolution rotations=solveRotations(x1,x2,c1,c2,maxAbsRotation);
      if (!rotations.valid) { valid=false; break; }
      c1=rotateZ(c1,rotations.dphi1); c2=rotateZ(c2,rotations.dphi2);
      totalR1+=rotations.dphi1; totalR2+=rotations.dphi2;
      done=iter+1;
    }

    outMassBefore=massPipi(p1,p2); outDiraBefore=dir;
    outDcaBefore=lineLineDca(x1,p1,x2,p2); outPt1=dpt1; outPt2=dpt2; outOpening=angle;
    outValid=valid?1:0; outIterations=done;

    if (valid)
    {
      outMassAfter=massPipi(c1,c2); outDiraAfter=pointing(add(c1,c2),flight);
      outDcaAfter=lineLineDca(x1,c1,x2,c2); outScale1=totalS1; outScale2=totalS2;
      outRot1=totalR1; outRot2=totalR2;
      hMassBefore.Fill(outMassBefore); hMassAfter.Fill(outMassAfter);
      hDiraBefore.Fill(outDiraBefore); hDiraAfter.Fill(outDiraAfter);
      hDcaBefore.Fill(outDcaBefore); hDcaAfter.Fill(outDcaAfter);
      hS1.Fill(outScale1); hS2.Fill(outScale2); hR1.Fill(outRot1); hR2.Fill(outRot2);
      hScales.Fill(outScale1,outScale2); hRotations.Fill(outRot1,outRot2);
      ++solved;
    }
    else
    {
      const float nan=std::numeric_limits<float>::quiet_NaN();
      outMassAfter=outDiraAfter=outDcaAfter=outScale1=outScale2=outRot1=outRot2=nan;
    }
    result.Fill();
  }

  output->Write();
  output->Close();
  std::cout<<"Selected: "<<selected<<" solved: "<<solved<<" wrote: "<<outputName<<std::endl;
}
