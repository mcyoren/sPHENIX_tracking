// Detailed QA of CM distortion maps for one "study" (a set of configurations).
// Writes a multi-page pdf; everything in mm (dR, R*dphi), empty map bins skipped.
//
//  1  title page: inputs and summary numbers
//  1b 1D projections in the layout of Ben Kimelman's slide (dphi [rad], dR [cm])
//  2  2D maps per configuration (dR / R dphi, South / North), common color scale
//  3  2D difference to the reference configuration (first in the list by default)
//  4  2D "higher order" residual: ring mean and first phi harmonic removed per ring
//  5  vs phi at six R rings, all configurations overlaid
//  6  vs R at eight phi values, all configurations overlaid
//  7  per-ring decomposition vs R: ring means <dR>, <R dphi> and first harmonic dx, dy
//  8  summary vs configuration: mean, RMS, RMS after harmonic removal, ring-averaged dx, dy
//
// inputs: *.list with lines "file label"
// root -l -b -q 'CMQA.C+("study.list","QA_study.pdf","title",0)'

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cmqa
{
  constexpr double kEmpty = -1e9;

  struct Ring
  {
    double R{0};
    int j{0};
    double A{0}, B{0}, dx{0}, dy{0};  // mm
    bool ok{false};
  };

  struct Config
  {
    std::string label;
    // [side][q]: q = 0 dR, q = 1 R dphi; mm; empty bins = kEmpty
    TH2D *h[2][2]{{nullptr, nullptr}, {nullptr, nullptr}};
    TH2D *hres[2][2]{{nullptr, nullptr}, {nullptr, nullptr}};  // higher-order residual
    std::vector<Ring> rings[2];
    double mean[2][2]{}, rms[2][2]{}, rmsH[2][2]{}, dxAvg[2]{}, dyAvg[2]{};
  };

  const char *qName[2] = {"#DeltaR [mm]", "R#Delta#phi [mm]"};
  const char *qShort[2] = {"dR", "Rdphi"};
  const char *sideName[2] = {"South", "North"};
  const char *sideExt[2] = {"_negz", "_posz"};

  bool filled(const TH2 *h, int i, int j) { return h->GetBinContent(i, j) > kEmpty / 2; }

  TH2D *makeMM(TH2 *hsrc, TH2 *hE, bool isPhi, const std::string &name)
  {
    auto *h = new TH2D(name.c_str(), "", hsrc->GetNbinsX(), hsrc->GetXaxis()->GetXmin(), hsrc->GetXaxis()->GetXmax(),
                       hsrc->GetNbinsY(), hsrc->GetYaxis()->GetXmin(), hsrc->GetYaxis()->GetXmax());
    h->SetDirectory(nullptr);
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
      for (int j = 1; j <= h->GetNbinsY(); ++j)
      {
        const bool guard = (i == 1 || j == 1 || i == h->GetNbinsX() || j == h->GetNbinsY());
        if (guard || hE->GetBinContent(i, j) <= 0)
        {
          h->SetBinContent(i, j, kEmpty);
          continue;
        }
        const double R = h->GetYaxis()->GetBinCenter(j);
        h->SetBinContent(i, j, isPhi ? 10 * R * hsrc->GetBinContent(i, j) : 10 * hsrc->GetBinContent(i, j));
      }
    }
    return h;
  }

  void analyse(Config &c)
  {
    for (int s = 0; s < 2; ++s)
    {
      TH2D *hR = c.h[s][0];
      TH2D *hP = c.h[s][1];
      for (int q = 0; q < 2; ++q)
      {
        c.hres[s][q] = static_cast<TH2D *>(c.h[s][q]->Clone(Form("%s_res", c.h[s][q]->GetName())));
        c.hres[s][q]->SetDirectory(nullptr);
      }
      double sum[2] = {0, 0}, sum2[2] = {0, 0}, sumH2[2] = {0, 0};
      int n = 0, nH = 0, nr = 0;
      for (int j = 2; j < hR->GetNbinsY(); ++j)
      {
        Ring ring;
        ring.j = j;
        ring.R = hR->GetYaxis()->GetBinCenter(j);
        std::vector<int> is;
        for (int i = 2; i < hR->GetNbinsX(); ++i)
        {
          if (!filled(hR, i, j))
          {
            continue;
          }
          is.push_back(i);
          for (int q = 0; q < 2; ++q)
          {
            const double v = c.h[s][q]->GetBinContent(i, j);
            sum[q] += v;
            sum2[q] += v * v;
          }
          ++n;
        }
        const int m = is.size();
        if (m < 12)
        {
          for (int i : is)
          {
            for (int q = 0; q < 2; ++q)
            {
              c.hres[s][q]->SetBinContent(i, j, kEmpty);
            }
          }
          c.rings[s].push_back(ring);
          continue;
        }
        TMatrixD M(2 * m, 4);
        TVectorD y(2 * m);
        for (int k = 0; k < m; ++k)
        {
          const double phi = hR->GetXaxis()->GetBinCenter(is[k]);
          const double co = std::cos(phi), si = std::sin(phi);
          M(k, 0) = 1; M(k, 1) = 0; M(k, 2) = co; M(k, 3) = si; y(k) = hR->GetBinContent(is[k], j);
          M(m + k, 0) = 0; M(m + k, 1) = 1; M(m + k, 2) = -si; M(m + k, 3) = co; y(m + k) = hP->GetBinContent(is[k], j);
        }
        TMatrixD Mt(TMatrixD::kTransposed, M);
        TMatrixD MtM = Mt * M;
        if (std::abs(MtM.Determinant()) < 1e-12)
        {
          c.rings[s].push_back(ring);
          continue;
        }
        TVectorD p = MtM.Invert() * (Mt * y);
        ring.A = p(0); ring.B = p(1); ring.dx = p(2); ring.dy = p(3); ring.ok = true;
        c.dxAvg[s] += p(2);
        c.dyAvg[s] += p(3);
        ++nr;
        for (int k = 0; k < m; ++k)
        {
          const double phi = hR->GetXaxis()->GetBinCenter(is[k]);
          const double co = std::cos(phi), si = std::sin(phi);
          const double rR = hR->GetBinContent(is[k], j) - p(0) - p(2) * co - p(3) * si;
          const double rP = hP->GetBinContent(is[k], j) - p(1) + p(2) * si - p(3) * co;
          c.hres[s][0]->SetBinContent(is[k], j, rR);
          c.hres[s][1]->SetBinContent(is[k], j, rP);
          // RMS after harmonic removal keeps the ring mean (space-charge-like part)
          sumH2[0] += (rR + p(0)) * (rR + p(0));
          sumH2[1] += (rP + p(1)) * (rP + p(1));
          ++nH;
        }
        c.rings[s].push_back(ring);
      }
      for (int q = 0; q < 2; ++q)
      {
        c.mean[s][q] = n ? sum[q] / n : 0;
        c.rms[s][q] = n ? std::sqrt(sum2[q] / n) : 0;
        c.rmsH[s][q] = nH ? std::sqrt(sumH2[q] / nH) : 0;
      }
      if (nr)
      {
        c.dxAvg[s] /= nr;
        c.dyAvg[s] /= nr;
      }
    }
  }

  int colorOf(int i, int n)
  {
    static const int fixed[] = {kBlack, kBlue, kRed, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet - 6, kGray + 2};
    if (n <= 9)
    {
      return fixed[i];
    }
    return TColor::GetColorPalette(i * (TColor::GetNumberOfColors() - 1) / (n - 1));
  }

  std::pair<int, int> grid(int n)
  {
    const int nx = static_cast<int>(std::ceil(std::sqrt(n * 1.5)));
    const int ny = (n + nx - 1) / nx;
    return {nx, ny};
  }

  void rangeOf(const std::vector<TH2D *> &hs, double &lo, double &hi, bool symmetric)
  {
    std::vector<double> v;
    for (auto *h : hs)
    {
      for (int i = 2; i < h->GetNbinsX(); ++i)
      {
        for (int j = 2; j < h->GetNbinsY(); ++j)
        {
          if (filled(h, i, j))
          {
            v.push_back(h->GetBinContent(i, j));
          }
        }
      }
    }
    if (v.empty())
    {
      lo = -1;
      hi = 1;
      return;
    }
    // 1% / 99% quantiles: single outlier bins do not set the scale
    std::sort(v.begin(), v.end());
    lo = v[static_cast<size_t>(0.01 * (v.size() - 1))];
    hi = v[static_cast<size_t>(0.99 * (v.size() - 1))];
    if (symmetric)
    {
      const double a = std::max(std::abs(lo), std::abs(hi));
      lo = -a;
      hi = a;
    }
    if (hi <= lo)
    {
      hi = lo + 1e-3;
    }
  }

  void draw2D(TH2D *h, const std::string &title, double lo, double hi)
  {
    h->SetTitle(Form("%s;#phi [rad];R [cm]", title.c_str()));
    h->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
    h->GetYaxis()->SetRangeUser(20, 80);
    h->SetMinimum(lo);
    h->SetMaximum(hi);
    h->SetTitleSize(0.06, "t");
    h->Draw("colz");
  }

  TGraph *profile(TH2D *h, bool alongPhi, int bin, int color, int istyle)
  {
    auto *g = new TGraph();
    if (alongPhi)
    {
      for (int i = 2; i < h->GetNbinsX(); ++i)
      {
        if (filled(h, i, bin))
        {
          g->AddPoint(h->GetXaxis()->GetBinCenter(i), h->GetBinContent(i, bin));
        }
      }
    }
    else
    {
      for (int j = 2; j < h->GetNbinsY(); ++j)
      {
        if (filled(h, bin, j))
        {
          g->AddPoint(h->GetYaxis()->GetBinCenter(j), h->GetBinContent(bin, j));
        }
      }
    }
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20 + istyle % 14);
    g->SetMarkerSize(0.6);
    g->SetLineWidth(2);
    return g;
  }

  // clears the canvas, writes a header line on top and returns the pad below it
  TPad *newPage(TCanvas *c, const std::string &header)
  {
    c->Clear();
    c->cd();
    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(62);
    lt.SetTextSize(0.022);
    lt.SetTextAlign(13);
    lt.DrawLatex(0.005, 0.995, header.c_str());
    auto *pg = new TPad(Form("pg_%p", static_cast<void *>(c)), "", 0, 0, 1, 0.965);
    pg->SetFillStyle(0);
    pg->Draw();
    pg->cd();
    return pg;
  }

  void drawMulti(TMultiGraph *mg, double x0, double x1)
  {
    if (!mg->GetListOfGraphs() || mg->GetListOfGraphs()->GetSize() == 0)
    {
      return;
    }
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(x0, x1);
    auto *l = new TLine(x0, 0, x1, 0);
    l->SetLineStyle(2);
    l->Draw();
  }
}  // namespace cmqa

void CMQA(const std::string &listfile, const std::string &outpdf, const std::string &title = "", int refIndex = 0)
{
  using namespace cmqa;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(64);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPaperSize(28, 20);  // landscape pdf page in cm (ROOT falls back to A4 portrait above ~A4)

  std::vector<Config> cfg;
  {
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
      ss >> file;
      std::getline(ss, label);
      label.erase(0, label.find_first_not_of(" \t"));
      TFile f(file.c_str());
      if (f.IsZombie())
      {
        std::cout << "cannot open " << file << std::endl;
        continue;
      }
      Config c;
      c.label = label.empty() ? file : label;
      bool ok = true;
      for (int s = 0; s < 2; ++s)
      {
        auto *hR = dynamic_cast<TH2 *>(f.Get(Form("hIntDistortionR%s", sideExt[s])));
        auto *hP = dynamic_cast<TH2 *>(f.Get(Form("hIntDistortionP%s", sideExt[s])));
        auto *hE = dynamic_cast<TH2 *>(f.Get(Form("hEntries%s", sideExt[s])));
        if (!hR || !hP || !hE)
        {
          ok = false;
          break;
        }
        c.h[s][0] = makeMM(hR, hE, false, Form("c%zu_s%d_dR", cfg.size(), s));
        c.h[s][1] = makeMM(hP, hE, true, Form("c%zu_s%d_Rdphi", cfg.size(), s));
      }
      if (!ok)
      {
        std::cout << "missing maps in " << file << std::endl;
        continue;
      }
      analyse(c);
      cfg.push_back(c);
    }
  }
  const int n = cfg.size();
  if (n == 0)
  {
    std::cout << "no inputs in " << listfile << std::endl;
    return;
  }
  refIndex = std::clamp(refIndex, 0, n - 1);
  gStyle->SetPalette(kRainBow);
  std::vector<int> col(n);
  for (int i = 0; i < n; ++i)
  {
    col[i] = colorOf(i, n);
  }
  gStyle->SetPalette(kBird);

  const std::string open = outpdf + "(";
  const std::string mid = outpdf;
  TCanvas *c = new TCanvas("cQA", "CM QA", 2400, 1400);

  // ---- 1: title + numbers ---------------------------------------------------
  {
    c->Clear();
    auto *pt = new TPaveText(0.02, 0.02, 0.98, 0.98, "NDC");
    pt->SetTextFont(82);
    pt->SetTextAlign(12);
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(std::min(0.025, 0.7 / (n + 12)));
    pt->AddText(Form("%s   (reference: %s)", title.c_str(), cfg[refIndex].label.c_str()));
    pt->AddText("all in mm, over filled map bins;  rmsH = RMS with the first phi harmonic removed per ring (ring mean kept);  dx, dy = ring-averaged first harmonic");
    pt->AddText(Form("%-30s | %7s %7s %7s %7s %7s %7s %6s %6s | %7s %7s %7s %7s %7s %7s %6s %6s", "configuration",
                     "S<dR>", "S rms", "S rmsH", "S<Rdp>", "S rms", "S rmsH", "S dx", "S dy",
                     "N<dR>", "N rms", "N rmsH", "N<Rdp>", "N rms", "N rmsH", "N dx", "N dy"));
    for (const auto &cf : cfg)
    {
      pt->AddText(Form("%-30s | %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %6.2f %6.2f | %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %6.2f %6.2f",
                       cf.label.c_str(),
                       cf.mean[0][0], cf.rms[0][0], cf.rmsH[0][0], cf.mean[0][1], cf.rms[0][1], cf.rmsH[0][1], cf.dxAvg[0], cf.dyAvg[0],
                       cf.mean[1][0], cf.rms[1][0], cf.rmsH[1][0], cf.mean[1][1], cf.rms[1][1], cf.rmsH[1][1], cf.dxAvg[1], cf.dyAvg[1]));
      printf("%-30s S: <dR> %7.3f rms %6.3f rmsH %6.3f <Rdp> %7.3f rms %6.3f rmsH %6.3f dx %6.2f dy %6.2f | N: <dR> %7.3f rms %6.3f rmsH %6.3f <Rdp> %7.3f rms %6.3f rmsH %6.3f dx %6.2f dy %6.2f\n",
             cf.label.c_str(),
             cf.mean[0][0], cf.rms[0][0], cf.rmsH[0][0], cf.mean[0][1], cf.rms[0][1], cf.rmsH[0][1], cf.dxAvg[0], cf.dyAvg[0],
             cf.mean[1][0], cf.rms[1][0], cf.rmsH[1][0], cf.mean[1][1], cf.rms[1][1], cf.rmsH[1][1], cf.dxAvg[1], cf.dyAvg[1]);
    }
    pt->Draw();
    c->Print(open.c_str(), "Title:summary");
  }

  // legend helper
  auto legendPad = [&](double textSize)
  {
    auto *leg = new TLegend(0.02, 0.02, 0.98, 0.98);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(textSize);
    if (n > 18)
    {
      leg->SetNColumns(2);
    }
    for (int i = 0; i < n; ++i)
    {
      auto *g = new TGraph();
      g->SetLineColor(col[i]);
      g->SetMarkerColor(col[i]);
      g->SetMarkerStyle(20 + i % 14);
      g->SetLineWidth(2);
      leg->AddEntry(g, cfg[i].label.c_str(), "lp");
    }
    leg->Draw();
  };

  // ---- 1b: Ben's "1D Distortion Projections" layout -----------------------
  //  dphi [rad] and dR [cm] vs phi at R = 35 cm, vs R at phi = 3pi/2; empty bins drawn as 0 as on his slide
  {
    TPad *pg = newPage(c, Form("%s  |  1D projections as on Ben's slide: vs #phi at R=35 cm, vs R at #phi=3#pi/2", title.c_str()));
    pg->Divide(5, 2, 0.002, 0.002);
    const char *ytit[2] = {"#DeltaR [cm]", "#Delta#phi [rad]"};
    for (int s = 0; s < 2; ++s)
    {
      for (int panel = 0; panel < 4; ++panel)
      {
        const bool vsPhi = panel < 2;
        const int q = (panel == 0 || panel == 2) ? 1 : 0;  // dphi first, then dR, as on the slide
        pg->cd(s * 5 + panel + 1);
        gPad->SetGridy();
        gPad->SetRightMargin(0.03);
        gPad->SetLeftMargin(0.16);
        std::vector<TH1D *> hs;
        double lo = 1e30, hi = -1e30;
        for (int i = 0; i < n; ++i)
        {
          const TH2D *h2 = cfg[i].h[s][q];
          const int fixedBin = vsPhi ? h2->GetYaxis()->FindBin(35.) : h2->GetXaxis()->FindBin(1.5 * TMath::Pi());
          const TAxis *ax = vsPhi ? h2->GetXaxis() : h2->GetYaxis();
          auto *h1 = new TH1D(Form("ben_%d_%d_%d", s, panel, i), "", ax->GetNbins(), ax->GetXmin(), ax->GetXmax());
          h1->SetDirectory(nullptr);
          for (int b = 1; b <= ax->GetNbins(); ++b)
          {
            const int ix = vsPhi ? b : fixedBin;
            const int iy = vsPhi ? fixedBin : b;
            double v = h2->GetBinContent(ix, iy);
            if (!filled(h2, ix, iy))
            {
              v = 0;
            }
            else
            {
              const double R = h2->GetYaxis()->GetBinCenter(iy);
              v = (q == 1) ? v / (10 * R) : v / 10;  // mm -> rad, mm -> cm
            }
            h1->SetBinContent(b, v);
            lo = std::min(lo, v);
            hi = std::max(hi, v);
          }
          h1->SetLineColor(col[i]);
          h1->SetLineWidth(2);
          h1->SetTitle(Form("%s vs %s, TPC %s;%s;%s", ytit[q], vsPhi ? "#phi" : "R", sideName[s], vsPhi ? "#phi [rad]" : "R [cm]", ytit[q]));
          h1->GetXaxis()->SetRangeUser(vsPhi ? 0 : 20, vsPhi ? 2 * TMath::Pi() : 80);
          hs.push_back(h1);
        }
        const double pad = 0.1 * std::max(hi - lo, 1e-6);
        for (size_t i = 0; i < hs.size(); ++i)
        {
          hs[i]->SetMinimum(lo - pad);
          hs[i]->SetMaximum(hi + pad);
          hs[i]->Draw(i == 0 ? "hist" : "hist same");
        }
      }
    }
    pg->cd(5);
    legendPad(0.06);
    c->Print(mid.c_str(), "Title:Ben-style projections");
  }

  // ---- 2-4: 2D pages -----------------------------------------------------
  auto [gx, gy] = grid(n);
  auto pages2D = [&](int kind)
  {
    // kind 0: maps, 1: difference to reference, 2: higher-order residual
    const char *kname[3] = {"map", Form("minus %s", cfg[refIndex].label.c_str()), "ring mean + 1st harmonic removed"};
    for (int s = 0; s < 2; ++s)
    {
      for (int q = 0; q < 2; ++q)
      {
        std::vector<TH2D *> hs;
        for (int i = 0; i < n; ++i)
        {
          TH2D *h = nullptr;
          if (kind == 0)
          {
            h = static_cast<TH2D *>(cfg[i].h[s][q]->Clone());
          }
          else if (kind == 2)
          {
            h = static_cast<TH2D *>(cfg[i].hres[s][q]->Clone());
          }
          else
          {
            if (i == refIndex)
            {
              continue;
            }
            h = static_cast<TH2D *>(cfg[i].h[s][q]->Clone());
            const TH2D *r = cfg[refIndex].h[s][q];
            for (int a = 1; a <= h->GetNbinsX(); ++a)
            {
              for (int b = 1; b <= h->GetNbinsY(); ++b)
              {
                h->SetBinContent(a, b, (filled(h, a, b) && filled(r, a, b)) ? h->GetBinContent(a, b) - r->GetBinContent(a, b) : kEmpty);
              }
            }
          }
          h->SetDirectory(nullptr);
          h->SetName(Form("p%d_%d_%d_%d", kind, s, q, i));
          h->SetTitle(cfg[i].label.c_str());
          hs.push_back(h);
        }
        if (hs.empty())
        {
          continue;
        }
        double lo, hi;
        rangeOf(hs, lo, hi, kind != 0);
        TPad *pg = newPage(c, Form("%s  |  %s, %s: %s", title.c_str(), qName[q], sideName[s], kname[kind]));
        auto [px, py] = grid(hs.size());
        pg->Divide(px, py, 0.002, 0.002);
        for (size_t k = 0; k < hs.size(); ++k)
        {
          pg->cd(k + 1);
          draw2D(hs[k], Form("%s", hs[k]->GetTitle()), lo, hi);
        }
        c->Print(mid.c_str(), Form("Title:2D %s %s %s", kind == 0 ? "map" : (kind == 1 ? "diff" : "residual"), qShort[q], sideName[s]));
      }
    }
  };
  (void) gx;
  (void) gy;
  pages2D(0);
  if (n > 1)
  {
    pages2D(1);
  }
  pages2D(2);

  // ---- 5: vs phi at six rings --------------------------------------------
  const double rings[6] = {33, 40, 48, 55, 63, 70};
  for (int q = 0; q < 2; ++q)
  {
    TPad *pg = newPage(c, Form("%s  |  %s vs #phi at fixed R", title.c_str(), qName[q]));
    pg->Divide(7, 2, 0.002, 0.002);
    for (int s = 0; s < 2; ++s)
    {
      for (int r = 0; r < 6; ++r)
      {
        pg->cd(s * 7 + r + 1);
        gPad->SetGridy();
        gPad->SetRightMargin(0.03);
        const int j = cfg[0].h[s][q]->GetYaxis()->FindBin(rings[r]);
        auto *mg = new TMultiGraph();
        mg->SetTitle(Form("%s, R=%.1f cm;#phi [rad];%s", sideName[s], cfg[0].h[s][q]->GetYaxis()->GetBinCenter(j), qName[q]));
        for (int i = 0; i < n; ++i)
        {
          auto *g = profile(cfg[i].h[s][q], true, j, col[i], i);
          if (g->GetN() > 0)
          {
            mg->Add(g, "LP");
          }
        }
        drawMulti(mg, 0, 2 * TMath::Pi());
      }
    }
    pg->cd(7);
    legendPad(0.06);
    c->Print(mid.c_str(), Form("Title:vs phi %s", qShort[q]));
  }

  // ---- 6: vs R at eight phi values ----------------------------------------
  for (int q = 0; q < 2; ++q)
  {
    TPad *pg = newPage(c, Form("%s  |  %s vs R at fixed #phi", title.c_str(), qName[q]));
    pg->Divide(5, 4, 0.002, 0.002);
    for (int s = 0; s < 2; ++s)
    {
      for (int k = 0; k < 8; ++k)
      {
        const double phi = (k + 0.5) * TMath::Pi() / 4;
        pg->cd(s * 10 + (k < 4 ? k : k + 1) + 1);
        gPad->SetGridy();
        gPad->SetRightMargin(0.03);
        const int i0 = cfg[0].h[s][q]->GetXaxis()->FindBin(phi);
        auto *mg = new TMultiGraph();
        mg->SetTitle(Form("%s, #phi=%.2f;R [cm];%s", sideName[s], cfg[0].h[s][q]->GetXaxis()->GetBinCenter(i0), qName[q]));
        for (int i = 0; i < n; ++i)
        {
          auto *g = profile(cfg[i].h[s][q], false, i0, col[i], i);
          if (g->GetN() > 0)
          {
            mg->Add(g, "LP");
          }
        }
        drawMulti(mg, 20, 80);
      }
    }
    pg->cd(5);
    legendPad(0.07);
    c->Print(mid.c_str(), Form("Title:vs R %s", qShort[q]));
  }

  // ---- 7: per-ring decomposition ------------------------------------------
  {
    TPad *pg = newPage(c, Form("%s  |  per-ring decomposition: ring means and first #phi harmonic", title.c_str()));
    pg->Divide(5, 2, 0.002, 0.002);
    const char *dname[4] = {"ring mean <#DeltaR> [mm]", "ring mean <R#Delta#phi> [mm]", "1st harmonic dx [mm]", "1st harmonic dy [mm]"};
    for (int s = 0; s < 2; ++s)
    {
      for (int d = 0; d < 4; ++d)
      {
        pg->cd(s * 5 + d + 1);
        gPad->SetGridy();
        gPad->SetRightMargin(0.03);
        auto *mg = new TMultiGraph();
        mg->SetTitle(Form("%s, %s;R [cm];%s", dname[d], sideName[s], dname[d]));
        for (int i = 0; i < n; ++i)
        {
          auto *g = new TGraph();
          for (const auto &ring : cfg[i].rings[s])
          {
            if (ring.ok)
            {
              const double v[4] = {ring.A, ring.B, ring.dx, ring.dy};
              g->AddPoint(ring.R, v[d]);
            }
          }
          g->SetLineColor(col[i]);
          g->SetMarkerColor(col[i]);
          g->SetMarkerStyle(20 + i % 14);
          g->SetMarkerSize(0.7);
          g->SetLineWidth(2);
          if (g->GetN() > 0)
          {
            mg->Add(g, "LP");
          }
        }
        drawMulti(mg, 20, 80);
      }
    }
    pg->cd(5);
    legendPad(0.06);
    c->Print(mid.c_str(), "Title:ring decomposition");
  }

  // ---- 8: summary vs configuration ----------------------------------------
  {
    TPad *pg = newPage(c, Form("%s  |  summary vs configuration (blue South, red North)", title.c_str()));
    pg->Divide(4, 2, 0.002, 0.002);
    const char *sname[8] = {"mean #DeltaR [mm]", "RMS #DeltaR [mm]", "mean R#Delta#phi [mm]", "RMS R#Delta#phi [mm]",
                            "first harmonic dx [mm]", "first harmonic dy [mm]", "RMS_{H} #DeltaR [mm]", "RMS_{H} R#Delta#phi [mm]"};
    for (int k = 0; k < 8; ++k)
    {
      pg->cd(k + 1);
      gPad->SetGridy();
      gPad->SetBottomMargin(0.35);
      gPad->SetRightMargin(0.03);
      TH1D *hs[2];
      double lo = 1e30, hi = -1e30;
      for (int s = 0; s < 2; ++s)
      {
        hs[s] = new TH1D(Form("sum_%d_%d", k, s), Form("%s;;%s", sname[k], sname[k]), n, 0, n);
        hs[s]->SetDirectory(nullptr);
        for (int i = 0; i < n; ++i)
        {
          const auto &cf = cfg[i];
          const double v[8] = {cf.mean[s][0], cf.rms[s][0], cf.mean[s][1], cf.rms[s][1], cf.dxAvg[s], cf.dyAvg[s], cf.rmsH[s][0], cf.rmsH[s][1]};
          hs[s]->SetBinContent(i + 1, v[k]);
          hs[s]->GetXaxis()->SetBinLabel(i + 1, cf.label.c_str());
          lo = std::min(lo, v[k]);
          hi = std::max(hi, v[k]);
        }
        hs[s]->SetMarkerStyle(s ? 24 : 20);
        hs[s]->SetMarkerSize(1.3);
        hs[s]->SetMarkerColor(s ? kRed : kBlue);
        hs[s]->SetLineColor(s ? kRed : kBlue);
        hs[s]->GetXaxis()->LabelsOption("v");
        hs[s]->GetXaxis()->SetLabelSize(n > 20 ? 0.035 : 0.05);
      }
      if (k == 1 || k == 3 || k >= 6)
      {
        lo = std::min(lo, 0.);
      }
      const double pad = 0.1 * std::max(hi - lo, 1e-3);
      hs[0]->SetMinimum(lo - pad);
      hs[0]->SetMaximum(hi + pad);
      hs[0]->Draw("P");
      hs[1]->Draw("P same");
      auto *l = new TLine(0, 0, n, 0);
      l->SetLineStyle(2);
      l->Draw();
      if (k == 0)
      {
        auto *leg = new TLegend(0.6, 0.8, 0.95, 0.93);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hs[0], "South", "p");
        leg->AddEntry(hs[1], "North", "p");
        leg->Draw();
      }
    }
    c->Print((outpdf + ")").c_str(), "Title:summary vs configuration");
  }
  std::cout << "wrote " << outpdf << std::endl;
}
