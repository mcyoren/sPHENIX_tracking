// CM distortion maps for several configurations, written to one multi-page pdf:
//   page 1: 1D projections as in Ben Kimelman's "1D Distortion Projections" slide,
//           dphi and dR vs phi at R = rProj and vs R at phi = phiProj, south/north
//   page 2: score per configuration (what is left after the Garfield correction):
//           mean and RMS of dR and R*dphi over all filled map bins, raw and after
//           removing the first phi harmonic per R ring (TPC tilt / offsets)
//   then  : one page of 2D maps (dR, dphi, south/north) per configuration,
//           same color scale for all configurations
//
// inputs: comma-separated "file:label" list, or a text file (*.list) with lines "file label"
// root -l -b -q 'PlotCMConfigProjections.C("maps.list","mapscan.pdf")'

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TStyle.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{
  std::vector<std::pair<std::string, std::string>> parseInputs(const std::string &inputs)
  {
    std::vector<std::pair<std::string, std::string>> out;
    if (inputs.size() > 5 && inputs.substr(inputs.size() - 5) == ".list")
    {
      std::ifstream in(inputs);
      std::string line;
      while (std::getline(in, line))
      {
        if (line.empty() || line[0] == '#')
        {
          continue;
        }
        std::istringstream ss(line);
        std::string file;
        std::string label;
        ss >> file;
        std::getline(ss, label);
        label.erase(0, label.find_first_not_of(" \t"));
        out.emplace_back(file, label.empty() ? file : label);
      }
      return out;
    }
    TObjArray *tokens = TString(inputs).Tokenize(",");
    for (int i = 0; i < tokens->GetEntries(); ++i)
    {
      TString item = static_cast<TObjString *>(tokens->At(i))->GetString();
      const int colon = item.Last(':');
      const TString file = colon > 0 ? TString(item(0, colon)) : item;
      const TString label = colon > 0 ? TString(item(colon + 1, item.Length())) : item;
      out.emplace_back(file.Data(), label.Data());
    }
    delete tokens;
    return out;
  }

  TH1D *projPhi(TH2 *h, double r, const std::string &name)
  {
    const int ib = h->GetYaxis()->FindBin(r);
    TH1D *p = h->ProjectionX(name.c_str(), ib, ib);
    p->SetDirectory(nullptr);
    p->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
    return p;
  }

  TH1D *projR(TH2 *h, double phi, const std::string &name)
  {
    const int ib = h->GetXaxis()->FindBin(phi);
    TH1D *p = h->ProjectionY(name.c_str(), ib, ib);
    p->SetDirectory(nullptr);
    p->GetXaxis()->SetRangeUser(20, 80);
    return p;
  }

  struct Score
  {
    double meanR{0};
    double rmsR{0};
    double meanP{0};
    double rmsP{0};
    double rmsRh{0};  // after removing the first harmonic per R ring
    double rmsPh{0};
    int n{0};
  };

  // dR in cm and dphi in rad -> scores in mm
  Score score(TH2 *hR, TH2 *hP, TH2 *hE)
  {
    Score s;
    double sumR = 0, sumR2 = 0, sumP = 0, sumP2 = 0, sumRh2 = 0, sumPh2 = 0;
    int nh = 0;
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
        const double dR = 10 * hR->GetBinContent(i, j);
        const double rdP = 10 * R * hP->GetBinContent(i, j);
        phis.push_back(hR->GetXaxis()->GetBinCenter(i));
        dRs.push_back(dR);
        rdPs.push_back(rdP);
        sumR += dR;
        sumR2 += dR * dR;
        sumP += rdP;
        sumP2 += rdP * rdP;
        ++s.n;
      }
      // first harmonic: dR = A + dx cos + dy sin, R dphi = B - dx sin + dy cos
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
        const double sn = std::sin(phis[k]);
        M(k, 0) = 1; M(k, 1) = 0; M(k, 2) = c; M(k, 3) = sn; y(k) = dRs[k];
        M(n + k, 0) = 0; M(n + k, 1) = 1; M(n + k, 2) = -sn; M(n + k, 3) = c; y(n + k) = rdPs[k];
      }
      TMatrixD Mt(TMatrixD::kTransposed, M);
      TMatrixD MtM = Mt * M;
      if (std::abs(MtM.Determinant()) < 1e-12)
      {
        continue;
      }
      TVectorD p = MtM.Invert() * (Mt * y);
      for (int k = 0; k < n; ++k)
      {
        // keep the ring mean (A, B): only the harmonic is removed
        const double c = std::cos(phis[k]);
        const double sn = std::sin(phis[k]);
        const double rR = dRs[k] - p(2) * c - p(3) * sn;
        const double rP = rdPs[k] + p(2) * sn - p(3) * c;
        sumRh2 += rR * rR;
        sumPh2 += rP * rP;
        ++nh;
      }
    }
    if (s.n > 0)
    {
      s.meanR = sumR / s.n;
      s.rmsR = std::sqrt(sumR2 / s.n);
      s.meanP = sumP / s.n;
      s.rmsP = std::sqrt(sumP2 / s.n);
    }
    if (nh > 0)
    {
      s.rmsRh = std::sqrt(sumRh2 / nh);
      s.rmsPh = std::sqrt(sumPh2 / nh);
    }
    return s;
  }
}  // namespace

void PlotCMConfigProjections(const std::string &inputs, const std::string &outfile = "CMConfigProjections.pdf",
                             double rProj = 35., double phiProj = 1.5 * TMath::Pi(), bool draw2D = true)
{
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(28, 20);  // landscape pdf page in cm (ROOT falls back to A4 portrait above ~A4)
  const auto files = parseInputs(inputs);
  const int nfiles = files.size();
  if (nfiles == 0)
  {
    std::cout << "no inputs" << std::endl;
    return;
  }

  // up to 9 inputs: distinct colors; more: palette gradient in input order
  const int fixedColors[] = {kGreen + 2, kBlue, kRed, kBlack, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet - 6, kGray + 2};
  gStyle->SetPalette(kRainBow);
  auto color = [&](int i)
  {
    if (nfiles <= 9)
    {
      return fixedColors[i];
    }
    return TColor::GetColorPalette(i * (TColor::GetNumberOfColors() - 1) / (nfiles - 1));
  };

  const char *vars[2] = {"P", "R"};
  const char *ytitles[2] = {"#Delta#phi [rad]", "#DeltaR [cm]"};
  const char *sides[2] = {"_negz", "_posz"};
  const char *sideNames[2] = {"TPC South", "TPC North"};

  // configs labelled "noGarfield" are references: drawn, but they do not set the axis ranges
  auto isRef = [](const std::string &label) { return label.find("noGarfield") != std::string::npos; };
  std::vector<bool> refOf;  // per projection entry index == file order

  std::vector<TH1D *> hists[2][4];
  std::vector<TH2 *> maps2D[2][2];  // [side][var]
  std::vector<Score> scores[2];
  std::vector<std::string> labels;

  TLegend *leg = new TLegend(0.05, 0.05, 0.95, 0.95);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(nfiles > 16 ? 2 : 1);

  for (int ifile = 0; ifile < nfiles; ++ifile)
  {
    const auto &[fname, label] = files[ifile];
    TFile *f = TFile::Open(fname.c_str());
    if (!f || f->IsZombie())
    {
      std::cout << "cannot open " << fname << std::endl;
      continue;
    }
    labels.push_back(label);
    refOf.push_back(isRef(label));

    for (int side = 0; side < 2; ++side)
    {
      auto *hE = dynamic_cast<TH2 *>(f->Get(Form("hEntries%s", sides[side])));
      TH2 *hv[2] = {nullptr, nullptr};
      for (int iv = 0; iv < 2; ++iv)
      {
        TH2 *h = dynamic_cast<TH2 *>(f->Get(Form("hIntDistortion%s%s", vars[iv], sides[side])));
        if (!h)
        {
          std::cout << "missing hIntDistortion" << vars[iv] << sides[side] << " in " << fname << std::endl;
          continue;
        }
        h->SetDirectory(nullptr);
        h->SetName(Form("%s_%d", h->GetName(), ifile));
        hv[iv] = h;
        maps2D[side][iv].push_back(h);
        for (int ip = 0; ip < 2; ++ip)
        {
          const std::string name = Form("h%s%s_%d_%s", vars[iv], sides[side], ifile, ip == 0 ? "phi" : "r");
          TH1D *p = (ip == 0) ? projPhi(h, rProj, name) : projR(h, phiProj, name);
          p->SetLineColor(color(ifile));
          p->SetLineWidth(2);
          p->SetTitle(Form("%s vs %s, %s", ytitles[iv], ip == 0 ? "#phi" : "R", sideNames[side]));
          p->GetXaxis()->SetTitle(ip == 0 ? "#phi [rad]" : "R [cm]");
          p->GetYaxis()->SetTitle(ytitles[iv]);
          hists[side][2 * iv + ip].push_back(p);
          if (side == 0 && iv == 0 && ip == 0)
          {
            leg->AddEntry(p, label.c_str(), "l");
          }
        }
      }
      scores[side].push_back((hv[0] && hv[1] && hE) ? score(hv[1], hv[0], hE) : Score{});
    }
    f->Close();
  }

  // ---- page 1: projections ------------------------------------------------
  TCanvas *c = new TCanvas("cCMConfigs", "CM configs", 2400, 800);
  c->Divide(5, 2);
  // columns: dphi(phi) dR(phi) dphi(R) dR(R) legend; slot = 2*iv+ip
  const int padOfSlot[4] = {1, 3, 2, 4};
  for (int side = 0; side < 2; ++side)
  {
    for (int slot = 0; slot < 4; ++slot)
    {
      auto &hv = hists[side][slot];
      if (hv.empty())
      {
        continue;
      }
      double ymin = 1e30;
      double ymax = -1e30;
      for (size_t i = 0; i < hv.size(); ++i)
      {
        if (refOf[i] && hv.size() > 1)
        {
          continue;
        }
        ymin = std::min(ymin, hv[i]->GetMinimum());
        ymax = std::max(ymax, hv[i]->GetMaximum());
      }
      const double pad = 0.1 * std::max(ymax - ymin, 1e-6);
      c->cd(side * 5 + padOfSlot[slot]);
      gPad->SetGridy();
      for (size_t i = 0; i < hv.size(); ++i)
      {
        hv[i]->SetMinimum(ymin - pad);
        hv[i]->SetMaximum(ymax + pad);
        hv[i]->Draw(i == 0 ? "hist" : "hist same");
      }
    }
  }
  c->cd(5);
  leg->Draw();
  c->SaveAs((outfile + "(").c_str());

  // ---- page 2: scores ---------------------------------------------------
  const int n = labels.size();
  printf("\n%-24s %-6s %8s %8s %8s %8s %9s %9s   (mm, filled map bins)\n",
         "config", "side", "<dR>", "rmsdR", "<Rdphi>", "rmsRdphi", "rmsdR_h1", "rmsRdp_h1");
  for (int i = 0; i < n; ++i)
  {
    for (int side = 0; side < 2; ++side)
    {
      const auto &s = scores[side][i];
      printf("%-24s %-6s %8.3f %8.3f %8.3f %8.3f %9.3f %9.3f\n", labels[i].c_str(), side ? "North" : "South",
             s.meanR, s.rmsR, s.meanP, s.rmsP, s.rmsRh, s.rmsPh);
    }
  }
  printf("  rms*_h1 = RMS after removing the first phi harmonic per R ring (tilt/offsets), ring mean kept\n\n");

  TCanvas *cs = new TCanvas("cScores", "scores", 2400, 900);
  cs->Divide(4, 2);
  const char *stitles[4] = {"RMS #DeltaR [mm]", "RMS R#Delta#phi [mm]", "mean #DeltaR [mm]", "mean R#Delta#phi [mm]"};
  for (int side = 0; side < 2; ++side)
  {
    for (int q = 0; q < 4; ++q)
    {
      cs->cd(side * 4 + q + 1);
      gPad->SetGridy();
      gPad->SetBottomMargin(0.3);
      auto *h = new TH1D(Form("hs_%d_%d", side, q), Form("%s, %s", stitles[q], sideNames[side]), n, 0, n);
      auto *hh = (q < 2) ? new TH1D(Form("hsh_%d_%d", side, q), "", n, 0, n) : nullptr;
      for (int i = 0; i < n; ++i)
      {
        const auto &s = scores[side][i];
        const double v[4] = {s.rmsR, s.rmsP, s.meanR, s.meanP};
        h->SetBinContent(i + 1, v[q]);
        h->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
        if (hh)
        {
          hh->SetBinContent(i + 1, q == 0 ? s.rmsRh : s.rmsPh);
        }
      }
      h->GetXaxis()->LabelsOption("v");
      h->GetXaxis()->SetLabelSize(n > 16 ? 0.03 : 0.045);
      h->SetMarkerStyle(20);
      h->SetMarkerColor(kBlue);
      h->SetLineColor(kBlue);
      double lo = 1e30;
      double hi = -1e30;
      for (int i = 0; i < n; ++i)
      {
        if (refOf[i] && n > 1)
        {
          continue;  // reference is off scale, read it from the printed table
        }
        lo = std::min(lo, h->GetBinContent(i + 1));
        hi = std::max(hi, h->GetBinContent(i + 1));
        if (hh)
        {
          lo = std::min(lo, hh->GetBinContent(i + 1));
          hi = std::max(hi, hh->GetBinContent(i + 1));
        }
      }
      if (q < 2)
      {
        lo = 0;
      }
      if (hh)
      {
        hh->SetMarkerStyle(24);
        hh->SetMarkerColor(kRed);
        hh->SetLineColor(kRed);
      }
      h->SetMinimum(lo - 0.1 * std::abs(hi - lo));
      h->SetMaximum(hi + 0.1 * std::abs(hi - lo));
      h->Draw("P");
      if (hh)
      {
        hh->Draw("P same");
        if (side == 0 && q == 0)
        {
          auto *l = new TLegend(0.15, 0.32, 0.85, 0.45);
          l->SetBorderSize(0);
          l->SetFillStyle(0);
          l->AddEntry(h, "raw", "p");
          l->AddEntry(hh, "1st harmonic removed", "p");
          l->Draw();
        }
      }
    }
  }
  cs->SaveAs((outfile + (draw2D ? "" : ")")).c_str());

  // ---- 2D maps, one page per configuration, common color scale -----------
  if (draw2D)
  {
    double zlo[2][2];
    double zhi[2][2];
    for (int side = 0; side < 2; ++side)
    {
      for (int iv = 0; iv < 2; ++iv)
      {
        zlo[side][iv] = 1e30;
        zhi[side][iv] = -1e30;
        for (size_t i = 0; i < maps2D[side][iv].size(); ++i)
        {
          if (refOf[i] && maps2D[side][iv].size() > 1)
          {
            continue;
          }
          TH2 *h = maps2D[side][iv][i];
          for (int ix = 2; ix < h->GetNbinsX(); ++ix)
          {
            for (int iy = 2; iy < h->GetNbinsY(); ++iy)
            {
              if (h->GetBinContent(ix, iy) != 0)
              {
                zlo[side][iv] = std::min(zlo[side][iv], h->GetBinContent(ix, iy));
                zhi[side][iv] = std::max(zhi[side][iv], h->GetBinContent(ix, iy));
              }
            }
          }
        }
      }
    }
    gStyle->SetPalette(kBird);
    TCanvas *c2 = new TCanvas("c2D", "2D maps", 1600, 1200);
    for (int i = 0; i < static_cast<int>(maps2D[0][0].size()); ++i)
    {
      c2->Clear();
      c2->Divide(2, 2);
      for (int iv = 1; iv >= 0; --iv)  // dR row first
      {
        for (int side = 0; side < 2; ++side)
        {
          c2->cd((1 - iv) * 2 + side + 1);
          gPad->SetRightMargin(0.15);
          TH2 *h = maps2D[side][iv][i];
          h->SetTitle(Form("%s: %s, %s;#phi [rad];R [cm]", labels[i].c_str(), ytitles[iv], sideNames[side]));
          h->GetXaxis()->SetRangeUser(0, 2 * TMath::Pi());
          h->GetYaxis()->SetRangeUser(20, 80);
          // empty bins below the color range are not painted
          for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
          {
            for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
            {
              if (h->GetBinContent(ix, iy) == 0)
              {
                h->SetBinContent(ix, iy, -1e9);
              }
            }
          }
          if (refOf[i])
          {
            h->SetMinimum(-1e8);  // own range: recomputed from filled bins below
            double lo = 1e30;
            double hi = -1e30;
            for (int ix = 2; ix < h->GetNbinsX(); ++ix)
            {
              for (int iy = 2; iy < h->GetNbinsY(); ++iy)
              {
                const double v = h->GetBinContent(ix, iy);
                if (v > -1e8)
                {
                  lo = std::min(lo, v);
                  hi = std::max(hi, v);
                }
              }
            }
            h->SetMinimum(lo);
            h->SetMaximum(hi);
          }
          else
          {
            h->SetMinimum(zlo[side][iv]);
            h->SetMaximum(zhi[side][iv]);
          }
          h->Draw("colz");
        }
      }
      const bool last = (i == static_cast<int>(maps2D[0][0].size()) - 1);
      c2->SaveAs((outfile + (last ? ")" : "")).c_str());
    }
  }
}
