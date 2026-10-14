// Survey-geometry check from CM distortion maps.
//
// A transverse offset d = (dx, dy) of the reconstructed CM pattern w.r.t. the
// truth pattern shows up as a first harmonic in phi:
//   dR(phi)       = A + dx cos(phi) + dy sin(phi)
//   R*dphi(phi)   = B - dx sin(phi) + dy cos(phi)
// For every R bin of hIntDistortion{R,P}_{negz,posz} this fits (A, B, dx, dy)
// over phi, and plots dx(R), dy(R) and A(R), B(R) per side for each input map.
//
// What to expect:
//  - a rigid TPC shift moves pads and CM together: invisible in the TPC frame
//  - a TPC tilt w.r.t. B is visible: electrons follow B, so the two sides get
//    opposite offsets of about (omega tau)^2/(1+(omega tau)^2) * tilt * 102 cm,
//    i.e. ~1.3 mm in x for rot_y = 1.48 mrad, ~0.3 mm in y for rot_x = 0.3 mrad
//  - comparing global-frame reco with the (TPC-frame) truth pattern adds the
//    survey transform at the readout plane: South ~(-3.2, -0.6) mm, North ~0
//
// root -l -b -q 'PlotCMSurveyCheck.C("a.root:label a,b.root:label b","survey_check.pdf")'

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TStyle.h>
#include <TVectorD.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  struct HarmonicFit
  {
    double R{0};
    double A{0};
    double B{0};
    double dx{0};
    double dy{0};
    int n{0};
  };

  // returns fits for all R bins with >= minPhiBins filled phi bins
  std::vector<HarmonicFit> fitSide(TH2 *hR, TH2 *hP, TH2 *hE, int minPhiBins = 12)
  {
    std::vector<HarmonicFit> fits;
    const int nphi = hR->GetNbinsX();
    for (int j = 2; j <= hR->GetNbinsY() - 1; ++j)
    {
      const double R = hR->GetYaxis()->GetBinCenter(j);
      std::vector<double> phis;
      std::vector<double> dRs;
      std::vector<double> rdPs;
      for (int i = 2; i <= nphi - 1; ++i)  // skip guard bins
      {
        if (hE->GetBinContent(i, j) <= 0)
        {
          continue;
        }
        phis.push_back(hR->GetXaxis()->GetBinCenter(i));
        dRs.push_back(hR->GetBinContent(i, j));
        rdPs.push_back(R * hP->GetBinContent(i, j));
      }
      if (static_cast<int>(phis.size()) < minPhiBins)
      {
        continue;
      }

      // joint linear least squares, parameters (A, B, dx, dy)
      const int n = phis.size();
      TMatrixD M(2 * n, 4);
      TVectorD y(2 * n);
      for (int k = 0; k < n; ++k)
      {
        const double c = std::cos(phis[k]);
        const double s = std::sin(phis[k]);
        M(k, 0) = 1;
        M(k, 1) = 0;
        M(k, 2) = c;
        M(k, 3) = s;
        y(k) = dRs[k];
        M(n + k, 0) = 0;
        M(n + k, 1) = 1;
        M(n + k, 2) = -s;
        M(n + k, 3) = c;
        y(n + k) = rdPs[k];
      }
      TMatrixD Mt(TMatrixD::kTransposed, M);
      TMatrixD MtM = Mt * M;
      if (std::abs(MtM.Determinant()) < 1e-12)
      {
        continue;
      }
      TVectorD p = MtM.Invert() * (Mt * y);
      fits.push_back({R, p(0), p(1), p(2), p(3), n});
    }
    return fits;
  }
}  // namespace

void PlotCMSurveyCheck(const std::string &inputs, const std::string &outfile = "CMSurveyCheck.pdf")
{
  const int colors[] = {kGreen + 2, kBlue, kRed, kBlack, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet - 6, kGray + 2};
  const int markers[] = {20, 21, 22, 23, 33, 34, 29, 47, 43};
  const int ncolors = sizeof(colors) / sizeof(colors[0]);
  const char *sides[2] = {"_negz", "_posz"};
  const char *sideNames[2] = {"South", "North"};

  // quantities: dx, dy [mm], A = <dR> [mm], B = <R dphi> [mm]
  const char *qnames[4] = {"dx [mm]", "dy [mm]", "#LT#DeltaR#GT [mm]", "#LTR#Delta#phi#GT [mm]"};
  TMultiGraph *mg[2][4];
  for (int s = 0; s < 2; ++s)
  {
    for (int q = 0; q < 4; ++q)
    {
      mg[s][q] = new TMultiGraph();
      mg[s][q]->SetTitle(Form("%s, %s;R [cm];%s", qnames[q], sideNames[s], qnames[q]));
    }
  }

  TLegend *leg = new TLegend(0.12, 0.3, 0.88, 0.75);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  printf("\n%-28s %-6s %9s %9s %9s %9s   (R-averaged, mm)\n", "input", "side", "dx", "dy", "<dR>", "<Rdphi>");

  TObjArray *tokens = TString(inputs).Tokenize(",");
  for (int ifile = 0; ifile < tokens->GetEntries(); ++ifile)
  {
    TString item = static_cast<TObjString *>(tokens->At(ifile))->GetString();
    const int colon = item.Last(':');
    const TString fname = colon > 0 ? TString(item(0, colon)) : item;
    const TString label = colon > 0 ? TString(item(colon + 1, item.Length())) : item;

    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie())
    {
      std::cout << "cannot open " << fname << std::endl;
      continue;
    }

    for (int s = 0; s < 2; ++s)
    {
      auto *hR = dynamic_cast<TH2 *>(f->Get(Form("hIntDistortionR%s", sides[s])));
      auto *hP = dynamic_cast<TH2 *>(f->Get(Form("hIntDistortionP%s", sides[s])));
      auto *hE = dynamic_cast<TH2 *>(f->Get(Form("hEntries%s", sides[s])));
      if (!hR || !hP || !hE)
      {
        std::cout << "missing maps for side " << sideNames[s] << " in " << fname << std::endl;
        continue;
      }
      const auto fits = fitSide(hR, hP, hE);
      if (fits.empty())
      {
        continue;
      }

      double sum[4] = {0, 0, 0, 0};
      TGraph *g[4];
      for (int q = 0; q < 4; ++q)
      {
        g[q] = new TGraph();
        g[q]->SetLineColor(colors[ifile % ncolors]);
        g[q]->SetMarkerColor(colors[ifile % ncolors]);
        g[q]->SetMarkerStyle(markers[ifile % ncolors]);
        g[q]->SetLineWidth(2);
      }
      for (const auto &fit : fits)
      {
        const double v[4] = {10 * fit.dx, 10 * fit.dy, 10 * fit.A, 10 * fit.B};  // cm -> mm
        for (int q = 0; q < 4; ++q)
        {
          g[q]->AddPoint(fit.R, v[q]);
          sum[q] += v[q];
        }
      }
      for (int q = 0; q < 4; ++q)
      {
        mg[s][q]->Add(g[q], "LP");
      }
      if (s == 0)
      {
        leg->AddEntry(g[0], label, "lp");
      }
      const double nf = fits.size();
      printf("%-28s %-6s %9.3f %9.3f %9.3f %9.3f\n", label.Data(), sideNames[s],
             sum[0] / nf, sum[1] / nf, sum[2] / nf, sum[3] / nf);
    }
    f->Close();
  }
  printf("\n");

  gStyle->SetPaperSize(28, 20);  // landscape pdf page in cm (ROOT falls back to A4 portrait above ~A4)
  TCanvas *c = new TCanvas("cCMSurvey", "CM survey check", 2000, 1000);
  c->Divide(4, 2);
  for (int s = 0; s < 2; ++s)
  {
    for (int q = 0; q < 4; ++q)
    {
      c->cd(s * 4 + q + 1);
      gPad->SetGridy();
      mg[s][q]->Draw("A");
      mg[s][q]->GetXaxis()->SetLimits(20, 80);
      TLine *l = new TLine(20, 0, 80, 0);
      l->SetLineStyle(2);
      l->Draw();
    }
  }
  c->cd(3);  // <dR> South: curves sit at the top and bottom of this pad
  leg->Draw();
  c->SaveAs(outfile.c_str());
}
