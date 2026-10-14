// Fit keff (relative to the CDB value) per space-charge map and side from a keff scan.
//
// Per R ring the map residual is split into ring means (A = <dR>, B = <R dphi>) and
// a first phi harmonic (tilt / offsets, not space charge) which is dropped. The ring
// means are linear in keff, so per ring  A_R(k) = a_R + s_R (k-1),  B_R(k) = b_R + t_R (k-1)
// from the scan points. Then keff is fitted by least squares over the rings:
//   M0: k only                     A_R(k), B_R(k) -> 0
//   M1: k + nuisances              A_R(k) + c, B_R(k) + phi0 R -> 0
//       c    = constant radial offset (CM pattern radius / readout radial offset)
//       phi0 = rigid rotation of the CM pattern relative to the pads
// Space charge has a definite R shape, the nuisances do not depend on keff, so M1
// lets the R shape (not the overall offsets) decide keff and the map.
//
// inputs: *.list with lines "file label", labels "<map>_k<factor with p>", e.g. alpha1p6_k0p975
// root -l -b -q 'FitCMKeff.C+("all.list","keff_fit.pdf")'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  // ring means (mm) after removing the first harmonic, keyed by R [cm] (rounded)
  std::map<int, std::pair<double, double>> ringMeans(TH2 *hR, TH2 *hP, TH2 *hE)
  {
    std::map<int, std::pair<double, double>> out;
    for (int j = 2; j <= hR->GetNbinsY() - 1; ++j)
    {
      const double R = hR->GetYaxis()->GetBinCenter(j);
      std::vector<double> phis, dRs, rdPs;
      for (int i = 2; i <= hR->GetNbinsX() - 1; ++i)
      {
        if (hE->GetBinContent(i, j) <= 0)
        {
          continue;
        }
        phis.push_back(hR->GetXaxis()->GetBinCenter(i));
        dRs.push_back(10 * hR->GetBinContent(i, j));
        rdPs.push_back(10 * R * hP->GetBinContent(i, j));
      }
      const int n = phis.size();
      if (n < 12)
      {
        continue;
      }
      TMatrixD M(2 * n, 4);
      TVectorD y(2 * n);
      for (int k = 0; k < n; ++k)
      {
        const double c = std::cos(phis[k]);
        const double s = std::sin(phis[k]);
        M(k, 0) = 1; M(k, 1) = 0; M(k, 2) = c; M(k, 3) = s; y(k) = dRs[k];
        M(n + k, 0) = 0; M(n + k, 1) = 1; M(n + k, 2) = -s; M(n + k, 3) = c; y(n + k) = rdPs[k];
      }
      TMatrixD Mt(TMatrixD::kTransposed, M);
      TMatrixD MtM = Mt * M;
      if (std::abs(MtM.Determinant()) < 1e-12)
      {
        continue;
      }
      TVectorD p = MtM.Invert() * (Mt * y);
      out[static_cast<int>(std::lround(R * 10))] = {p(0), p(1)};
    }
    return out;
  }

  struct FitResult
  {
    double k{std::nan("")};
    double c{0};
    double phi0{0};  // mrad
    double rms{0};   // mm, over A and B of all rings
    double rmsA{0};
    double rmsB{0};
    int nring{0};
  };

  // rows: per ring two equations; unknowns u = k-1 (+ c, phi0 if nuisance)
  FitResult fitK(const std::vector<double> &R, const std::vector<double> &a, const std::vector<double> &s,
                 const std::vector<double> &b, const std::vector<double> &t, bool nuisance)
  {
    FitResult res;
    const int n = R.size();
    const int np = nuisance ? 3 : 1;
    if (n < np + 1)
    {
      return res;
    }
    TMatrixD M(2 * n, np);
    TVectorD y(2 * n);
    for (int i = 0; i < n; ++i)
    {
      M(i, 0) = s[i];
      M(n + i, 0) = t[i];
      if (nuisance)
      {
        M(i, 1) = 1;
        M(i, 2) = 0;
        M(n + i, 1) = 0;
        M(n + i, 2) = 0.01 * R[i];  // phi0 [mrad] * R [cm] -> mm
      }
      y(i) = -a[i];
      y(n + i) = -b[i];
    }
    TMatrixD Mt(TMatrixD::kTransposed, M);
    TMatrixD MtM = Mt * M;
    if (std::abs(MtM.Determinant()) < 1e-15)
    {
      return res;
    }
    TVectorD p = MtM.Invert() * (Mt * y);
    res.k = 1 + p(0);
    if (nuisance)
    {
      res.c = p(1);
      res.phi0 = p(2);
    }
    double sA = 0, sB = 0;
    for (int i = 0; i < n; ++i)
    {
      const double rA = a[i] + s[i] * p(0) + (nuisance ? p(1) : 0);
      const double rB = b[i] + t[i] * p(0) + (nuisance ? 0.01 * R[i] * p(2) : 0);
      sA += rA * rA;
      sB += rB * rB;
    }
    res.rmsA = std::sqrt(sA / n);
    res.rmsB = std::sqrt(sB / n);
    res.rms = std::sqrt((sA + sB) / (2 * n));
    res.nring = n;
    return res;
  }
}  // namespace

void FitCMKeff(const std::string &listfile, const std::string &outfile = "keff_fit.pdf")
{
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(28, 20);  // landscape pdf page in cm (ROOT falls back to A4 portrait above ~A4)
  const char *sides[2] = {"_negz", "_posz"};
  const char *sideNames[2] = {"South", "North"};

  // map name -> k -> ring means [side]
  std::map<std::string, std::map<double, std::map<int, std::pair<double, double>>[2]>> data;
  std::vector<std::string> mapOrder;

  std::ifstream in(listfile);
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    std::istringstream ss(line);
    std::string file, label;
    ss >> file >> label;
    const auto pos = label.rfind("_k");
    if (pos == std::string::npos)
    {
      continue;  // e.g. noGarfield
    }
    const std::string mapname = label.substr(0, pos);
    std::string kstr = label.substr(pos + 2);
    std::replace(kstr.begin(), kstr.end(), 'p', '.');
    const double k = std::stod(kstr);

    TFile f(file.c_str());
    if (f.IsZombie())
    {
      std::cout << "cannot open " << file << std::endl;
      continue;
    }
    if (std::find(mapOrder.begin(), mapOrder.end(), mapname) == mapOrder.end())
    {
      mapOrder.push_back(mapname);
    }
    for (int side = 0; side < 2; ++side)
    {
      auto *hR = dynamic_cast<TH2 *>(f.Get(Form("hIntDistortionR%s", sides[side])));
      auto *hP = dynamic_cast<TH2 *>(f.Get(Form("hIntDistortionP%s", sides[side])));
      auto *hE = dynamic_cast<TH2 *>(f.Get(Form("hEntries%s", sides[side])));
      if (hR && hP && hE)
      {
        data[mapname][k][side] = ringMeans(hR, hP, hE);
      }
    }
  }

  printf("\n%-10s %-6s | %-22s | %-48s\n", "map", "side", "M0: k only", "M1: k + const dR + rigid rotation");
  printf("%-10s %-6s | %7s %13s | %7s %8s %10s %8s %8s %8s\n", "", "", "k/kCDB", "rms[mm]", "k/kCDB", "c[mm]", "phi0[mrad]",
         "rms[mm]", "rmsdR", "rmsRdphi");

  const int colors[] = {kBlack, kBlue, kGreen + 2, kRed, kMagenta + 1, kOrange + 7};
  TMultiGraph *mg[2][2];  // [side][A/B] residual ring means at the M1 best fit
  TMultiGraph *mg1[2][2]; // [side][A/B] ring means at k = 1, no nuisances
  for (int side = 0; side < 2; ++side)
  {
    for (int q = 0; q < 2; ++q)
    {
      mg[side][q] = new TMultiGraph();
      mg[side][q]->SetTitle(Form("%s after M1 fit, %s;R [cm];%s [mm]", q ? "<R#Delta#phi>" : "<#DeltaR>", sideNames[side], q ? "<R#Delta#phi>" : "<#DeltaR>"));
      mg1[side][q] = new TMultiGraph();
      mg1[side][q]->SetTitle(Form("%s at k = k_{CDB}, %s;R [cm];%s [mm]", q ? "<R#Delta#phi>" : "<#DeltaR>", sideNames[side], q ? "<R#Delta#phi>" : "<#DeltaR>"));
    }
  }
  TLegend *leg = new TLegend(0.15, 0.6, 0.6, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  for (size_t im = 0; im < mapOrder.size(); ++im)
  {
    const auto &mapname = mapOrder[im];
    const auto &scan = data[mapname];
    for (int side = 0; side < 2; ++side)
    {
      // rings present at every k
      std::map<int, int> count;
      for (const auto &[k, sd] : scan)
      {
        for (const auto &[r, v] : sd[side])
        {
          ++count[r];
        }
      }
      std::vector<double> R, a, s, b, t;
      for (const auto &[r, nk] : count)
      {
        if (nk != static_cast<int>(scan.size()) || nk < 2)
        {
          continue;
        }
        // linear regression in u = k-1
        double su = 0, suu = 0, sA = 0, suA = 0, sB = 0, suB = 0;
        for (const auto &[k, sd] : scan)
        {
          const double u = k - 1;
          const auto &v = sd[side].at(r);
          su += u; suu += u * u; sA += v.first; suA += u * v.first; sB += v.second; suB += u * v.second;
        }
        const double n = nk;
        const double den = n * suu - su * su;
        const double slA = (n * suA - su * sA) / den;
        const double slB = (n * suB - su * sB) / den;
        R.push_back(0.1 * r);
        a.push_back((sA - slA * su) / n);
        s.push_back(slA);
        b.push_back((sB - slB * su) / n);
        t.push_back(slB);
      }
      const FitResult f0 = fitK(R, a, s, b, t, false);
      const FitResult f1 = fitK(R, a, s, b, t, true);
      printf("%-10s %-6s | %7.3f %13.3f | %7.3f %8.3f %10.3f %8.3f %8.3f %8.3f\n", mapname.c_str(), sideNames[side],
             f0.k, f0.rms, f1.k, f1.c, f1.phi0, f1.rms, f1.rmsA, f1.rmsB);

      TGraph *g[2] = {new TGraph(), new TGraph()};
      TGraph *g1[2] = {new TGraph(), new TGraph()};
      for (size_t i = 0; i < R.size(); ++i)
      {
        g[0]->AddPoint(R[i], a[i] + s[i] * (f1.k - 1) + f1.c);
        g[1]->AddPoint(R[i], b[i] + t[i] * (f1.k - 1) + 0.01 * R[i] * f1.phi0);
        g1[0]->AddPoint(R[i], a[i]);
        g1[1]->AddPoint(R[i], b[i]);
      }
      for (int q = 0; q < 2; ++q)
      {
        for (TGraph *gg : {g[q], g1[q]})
        {
          gg->SetLineColor(colors[im % 6]);
          gg->SetMarkerColor(colors[im % 6]);
          gg->SetMarkerStyle(20 + im);
          gg->SetLineWidth(2);
        }
        mg[side][q]->Add(g[q], "LP");
        mg1[side][q]->Add(g1[q], "LP");
      }
      if (side == 0)
      {
        leg->AddEntry(g[0], mapname.c_str(), "lp");
      }
    }
  }
  printf("  k/kCDB: fitted keff relative to the CDB keff; rms: what the fit leaves in the ring means\n\n");

  TCanvas *c = new TCanvas("cKeff", "keff fit", 2000, 1000);
  c->Divide(4, 2);
  for (int side = 0; side < 2; ++side)
  {
    for (int q = 0; q < 2; ++q)
    {
      for (int which = 0; which < 2; ++which)
      {
        c->cd(side * 4 + which * 2 + q + 1);
        gPad->SetGridy();
        TMultiGraph *m = which ? mg[side][q] : mg1[side][q];
        m->Draw("A");
        auto *l = new TLine(m->GetXaxis()->GetXmin(), 0, m->GetXaxis()->GetXmax(), 0);
        l->SetLineStyle(2);
        l->Draw();
      }
    }
  }
  c->cd(1);
  leg->Draw();
  c->SaveAs(outfile.c_str());
}
