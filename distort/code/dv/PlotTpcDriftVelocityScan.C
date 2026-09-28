// PlotTpcDriftVelocityScan.C
//
// Post-process tpc_dv_scan_run*.root files produced by
// Fun4All_TpcDriftVelocityScan.C.
//
// Usage:
//   root --web=off -l -b -q 'PlotTpcDriftVelocityScan.C()'
//
// Optional:
//   root --web=off -l -b -q \
//     'PlotTpcDriftVelocityScan.C("/sphenix/tg/tg01/hf/mitrankov/dv/rootfiles",
//                                 "/sphenix/tg/tg01/hf/mitrankov/dv/plots")'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  struct RunData
  {
    int run = 0;
    double cdbDV_cm_ns = std::numeric_limits<double>::quiet_NaN();
    double cdbT0_ns = std::numeric_limits<double>::quiet_NaN();
    double correctedDV_cm_ns = std::numeric_limits<double>::quiet_NaN();
    double garfieldSouth_cm_us = std::numeric_limits<double>::quiet_NaN();
    double garfieldNorth_cm_us = std::numeric_limits<double>::quiet_NaN();
    double garfieldCombined_cm_us = std::numeric_limits<double>::quiet_NaN();
    double cmField_Vcm = std::numeric_limits<double>::quiet_NaN();
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperatureK = std::numeric_limits<double>::quiet_NaN();

    std::unique_ptr<TProfile> vzRSouth;
    std::unique_ptr<TProfile> vzRNorth;
    std::unique_ptr<TProfile> vzZSouth;
    std::unique_ptr<TProfile> vzZNorth;
    std::unique_ptr<TProfile> trajRSouth;
    std::unique_ptr<TProfile> trajRNorth;
    std::unique_ptr<TProfile2D> vzRZSouth;
    std::unique_ptr<TProfile2D> vzRZNorth;
    std::unique_ptr<TProfile2D> vzRPhiSouth;
    std::unique_ptr<TProfile2D> vzRPhiNorth;
  };

  template <class T>
  std::unique_ptr<T> CloneDetached(TFile *f, const char *name, const int run)
  {
    auto *obj = dynamic_cast<T *>(f->Get(name));
    if (!obj) return nullptr;

    auto *copy = dynamic_cast<T *>(obj->Clone(Form("%s_run%d_copy", name, run)));
    if (!copy) return nullptr;
    copy->SetDirectory(nullptr);
    return std::unique_ptr<T>(copy);
  }

  bool finitePositive(const double x)
  {
    return std::isfinite(x) && x > 0.0;
  }

  void SaveCanvas(TCanvas &c, const std::string &plotDir, const std::string &base)
  {
    c.SaveAs((plotDir + "/" + base + ".png").c_str());
    c.SaveAs((plotDir + "/" + base + ".pdf").c_str());
  }

  void StyleGraph(TGraph *g, int markerStyle)
  {
    if (!g) return;
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2);
  }

  std::unique_ptr<TGraph> MakeRunGraph(
      const std::vector<RunData> &runs,
      double RunData::*member,
      const char *name)
  {
    auto g = std::make_unique<TGraph>();
    g->SetName(name);
    for (const auto &r : runs)
    {
      const double y = r.*member;
      if (std::isfinite(y)) g->SetPoint(g->GetN(), r.run, y);
    }
    return g;
  }

  std::unique_ptr<TGraph> MakeRatioRunGraph(
      const std::vector<RunData> &runs,
      double RunData::*garfieldMember,
      const char *name)
  {
    auto g = std::make_unique<TGraph>();
    g->SetName(name);
    for (const auto &r : runs)
    {
      const double garf = r.*garfieldMember;
      const double cdb = 1000.0 * r.cdbDV_cm_ns;
      if (std::isfinite(garf) && finitePositive(cdb))
      {
        g->SetPoint(g->GetN(), r.run, garf / cdb);
      }
    }
    return g;
  }

  std::unique_ptr<TGraph> MakeRelativeDifferenceRunGraph(
      const std::vector<RunData> &runs,
      double RunData::*garfieldMember,
      const char *name)
  {
    auto g = std::make_unique<TGraph>();
    g->SetName(name);
    for (const auto &r : runs)
    {
      const double garf = r.*garfieldMember;
      const double cdb = 1000.0 * r.cdbDV_cm_ns;
      if (std::isfinite(garf) && finitePositive(cdb))
      {
        g->SetPoint(g->GetN(), r.run, 100.0 * (garf - cdb) / cdb);
      }
    }
    return g;
  }

  // Average profile/CDB ratio over all runs, with run-to-run RMS/sqrt(N) as error.
  std::unique_ptr<TGraphErrors> AverageProfileRatio(
      const std::vector<RunData> &runs,
      std::unique_ptr<TProfile> RunData::*profileMember,
      const char *name)
  {
    const TProfile *templ = nullptr;
    for (const auto &r : runs)
    {
      const auto &p = r.*profileMember;
      if (p)
      {
        templ = p.get();
        break;
      }
    }
    if (!templ) return nullptr;

    auto g = std::make_unique<TGraphErrors>();
    g->SetName(name);

    for (int ib = 1; ib <= templ->GetNbinsX(); ++ib)
    {
      std::vector<double> vals;
      double xsum = 0.0;
      int nx = 0;

      for (const auto &r : runs)
      {
        const auto &p = r.*profileMember;
        const double cdb = 1000.0 * r.cdbDV_cm_ns;
        if (!p || !finitePositive(cdb)) continue;
        if (p->GetBinEntries(ib) <= 0) continue;

        const double v = p->GetBinContent(ib);
        if (!std::isfinite(v)) continue;

        vals.push_back(v / cdb);
        xsum += p->GetXaxis()->GetBinCenter(ib);
        ++nx;
      }

      if (vals.empty()) continue;

      double mean = 0.0;
      for (double v : vals) mean += v;
      mean /= vals.size();

      double var = 0.0;
      for (double v : vals) var += (v - mean) * (v - mean);
      if (vals.size() > 1) var /= (vals.size() - 1);

      const double err = vals.size() > 1 ? std::sqrt(var / vals.size()) : 0.0;
      const double x = nx > 0 ? xsum / nx : templ->GetXaxis()->GetBinCenter(ib);

      const int ip = g->GetN();
      g->SetPoint(ip, x, mean);
      g->SetPointError(ip, 0.0, err);
    }

    return g;
  }

  std::unique_ptr<TH2D> AverageProfile2DRatio(
      const std::vector<RunData> &runs,
      std::unique_ptr<TProfile2D> RunData::*profileMember,
      const char *name,
      const char *title)
  {
    const TProfile2D *templ = nullptr;
    for (const auto &r : runs)
    {
      const auto &p = r.*profileMember;
      if (p)
      {
        templ = p.get();
        break;
      }
    }
    if (!templ) return nullptr;

    const int nx = templ->GetNbinsX();
    const int ny = templ->GetNbinsY();

    auto h = std::make_unique<TH2D>(
        name, title,
        nx, templ->GetXaxis()->GetXmin(), templ->GetXaxis()->GetXmax(),
        ny, templ->GetYaxis()->GetXmin(), templ->GetYaxis()->GetXmax());

    for (int ix = 1; ix <= nx; ++ix)
    {
      for (int iy = 1; iy <= ny; ++iy)
      {
        double sum = 0.0;
        int n = 0;

        for (const auto &r : runs)
        {
          const auto &p = r.*profileMember;
          const double cdb = 1000.0 * r.cdbDV_cm_ns;
          if (!p || !finitePositive(cdb)) continue;

          const int bin = p->GetBin(ix, iy);
          if (p->GetBinEntries(bin) <= 0) continue;

          const double v = p->GetBinContent(bin);
          if (!std::isfinite(v)) continue;

          sum += v / cdb;
          ++n;
        }

        if (n > 0) h->SetBinContent(ix, iy, sum / n);
      }
    }

    return h;
  }

  std::unique_ptr<TH2D> RunProfileRatioMap(
      const std::vector<RunData> &runs,
      std::unique_ptr<TProfile> RunData::*profileMember,
      const char *name,
      const char *title,
      const char *yTitle)
  {
    const TProfile *templ = nullptr;
    for (const auto &r : runs)
    {
      const auto &p = r.*profileMember;
      if (p)
      {
        templ = p.get();
        break;
      }
    }
    if (!templ || runs.empty()) return nullptr;

    const int nx = static_cast<int>(runs.size());
    const int ny = templ->GetNbinsX();

    auto h = std::make_unique<TH2D>(
        name, title,
        nx, 0.5, nx + 0.5,
        ny, templ->GetXaxis()->GetXmin(), templ->GetXaxis()->GetXmax());

    h->GetXaxis()->SetTitle("run");
    h->GetYaxis()->SetTitle(yTitle);
    h->GetZaxis()->SetTitle("<v_{z,out}> / v_{CDB}");

    for (int ir = 0; ir < nx; ++ir)
    {
      const auto &r = runs[ir];
      const auto &p = r.*profileMember;
      const double cdb = 1000.0 * r.cdbDV_cm_ns;

      h->GetXaxis()->SetBinLabel(ir + 1, Form("%d", r.run));

      if (!p || !finitePositive(cdb)) continue;

      for (int ib = 1; ib <= ny; ++ib)
      {
        if (p->GetBinEntries(ib) <= 0) continue;
        const double v = p->GetBinContent(ib);
        if (std::isfinite(v)) h->SetBinContent(ir + 1, ib, v / cdb);
      }
    }

    h->LabelsOption("v", "X");
    return h;
  }
}

void PlotTpcDriftVelocityScan(
    const std::string &rootDir = "/sphenix/tg/tg01/hf/mitrankov/dv/rootfiles",
    const std::string &plotDir = "/sphenix/tg/tg01/hf/mitrankov/dv/plots")
{
  gSystem->mkdir(plotDir.c_str(), true);
  gStyle->SetOptStat(0);

  std::vector<std::filesystem::path> files;
  for (const auto &entry : std::filesystem::directory_iterator(rootDir))
  {
    if (!entry.is_regular_file()) continue;
    const auto name = entry.path().filename().string();
    if (name.rfind("tpc_dv_scan_run", 0) == 0 && entry.path().extension() == ".root")
      files.push_back(entry.path());
  }
  std::sort(files.begin(), files.end());

  std::vector<RunData> runs;

  for (const auto &path : files)
  {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
      std::cerr << "Skipping unreadable file: " << path << std::endl;
      continue;
    }

    auto *tree = dynamic_cast<TTree *>(f->Get("run_summary"));
    if (!tree || tree->GetEntries() < 1)
    {
      std::cerr << "Skipping file without run_summary: " << path << std::endl;
      continue;
    }

    RunData d;

    tree->SetBranchAddress("run", &d.run);
    tree->SetBranchAddress("cdb_drift_velocity_cm_ns", &d.cdbDV_cm_ns);
    tree->SetBranchAddress("cdb_tzero_ns", &d.cdbT0_ns);
    tree->SetBranchAddress("cdb_drift_velocity_t0corr_cm_ns", &d.correctedDV_cm_ns);
    tree->SetBranchAddress("garfield_mean_vz_south_cm_us", &d.garfieldSouth_cm_us);
    tree->SetBranchAddress("garfield_mean_vz_north_cm_us", &d.garfieldNorth_cm_us);
    tree->SetBranchAddress("garfield_mean_vz_combined_cm_us", &d.garfieldCombined_cm_us);
    tree->SetBranchAddress("cm_field_Vcm", &d.cmField_Vcm);
    tree->SetBranchAddress("pressure", &d.pressure);
    tree->SetBranchAddress("temperature_K", &d.temperatureK);
    tree->GetEntry(0);

    d.vzRSouth = CloneDetached<TProfile>(f.get(), "p_vz_vs_r_south", d.run);
    d.vzRNorth = CloneDetached<TProfile>(f.get(), "p_vz_vs_r_north", d.run);
    d.vzZSouth = CloneDetached<TProfile>(f.get(), "p_vz_vs_z_south", d.run);
    d.vzZNorth = CloneDetached<TProfile>(f.get(), "p_vz_vs_z_north", d.run);
    d.trajRSouth = CloneDetached<TProfile>(f.get(), "p_trajectory_vz_vs_r_south", d.run);
    d.trajRNorth = CloneDetached<TProfile>(f.get(), "p_trajectory_vz_vs_r_north", d.run);
    d.vzRZSouth = CloneDetached<TProfile2D>(f.get(), "p_vz_rz_south", d.run);
    d.vzRZNorth = CloneDetached<TProfile2D>(f.get(), "p_vz_rz_north", d.run);
    d.vzRPhiSouth = CloneDetached<TProfile2D>(f.get(), "p_vz_rphi_south", d.run);
    d.vzRPhiNorth = CloneDetached<TProfile2D>(f.get(), "p_vz_rphi_north", d.run);

    runs.push_back(std::move(d));
  }

  std::sort(runs.begin(), runs.end(),
            [](const RunData &a, const RunData &b) { return a.run < b.run; });

  std::cout << "Loaded " << runs.size() << " run files from " << rootDir << std::endl;
  if (runs.empty()) return;

  // ----------------------------------------------------------------------
  // 1. Mean PHGarfield drift velocity vs CDB
  // ----------------------------------------------------------------------
  auto gSouth = MakeRunGraph(runs, &RunData::garfieldSouth_cm_us, "g_mean_south");
  auto gNorth = MakeRunGraph(runs, &RunData::garfieldNorth_cm_us, "g_mean_north");
  auto gCombined = MakeRunGraph(runs, &RunData::garfieldCombined_cm_us, "g_mean_combined");

  auto gCdb = std::make_unique<TGraph>();
  gCdb->SetName("g_cdb");
  auto gCdbT0 = std::make_unique<TGraph>();
  gCdbT0->SetName("g_cdb_t0corr");

  for (const auto &r : runs)
  {
    if (finitePositive(r.cdbDV_cm_ns))
      gCdb->SetPoint(gCdb->GetN(), r.run, 1000.0 * r.cdbDV_cm_ns);

    if (finitePositive(r.correctedDV_cm_ns))
      gCdbT0->SetPoint(gCdbT0->GetN(), r.run, 1000.0 * r.correctedDV_cm_ns);
  }

  StyleGraph(gSouth.get(), 20);
  StyleGraph(gNorth.get(), 21);
  StyleGraph(gCombined.get(), 24);
  StyleGraph(gCdb.get(), 25);
  StyleGraph(gCdbT0.get(), 26);

  {
    TCanvas c("c_mean_vs_run", "", 1600, 900);
    c.SetGrid();
    gCombined->SetTitle("Mean drift velocity by run;run;v_{z} [cm/#mus]");
    gCombined->Draw("APL");
    gSouth->Draw("PL SAME");
    gNorth->Draw("PL SAME");
    gCdb->Draw("PL SAME");
    if (gCdbT0->GetN() > 0) gCdbT0->Draw("PL SAME");

    TLegend leg(0.14, 0.69, 0.39, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(gCombined.get(), "PHGarfield combined", "lp");
    leg.AddEntry(gSouth.get(), "PHGarfield south", "lp");
    leg.AddEntry(gNorth.get(), "PHGarfield north", "lp");
    leg.AddEntry(gCdb.get(), "CDB drift velocity", "lp");
    if (gCdbT0->GetN() > 0) leg.AddEntry(gCdbT0.get(), "CDB t0-corrected", "lp");
    leg.Draw();

    SaveCanvas(c, plotDir, "01_mean_drift_velocity_vs_run");
  }

  // ----------------------------------------------------------------------
  // 2. PHGarfield/CDB ratio and relative difference vs run
  // ----------------------------------------------------------------------
  auto gRatioS = MakeRatioRunGraph(runs, &RunData::garfieldSouth_cm_us, "g_ratio_south");
  auto gRatioN = MakeRatioRunGraph(runs, &RunData::garfieldNorth_cm_us, "g_ratio_north");
  auto gRatioC = MakeRatioRunGraph(runs, &RunData::garfieldCombined_cm_us, "g_ratio_combined");
  StyleGraph(gRatioS.get(), 20);
  StyleGraph(gRatioN.get(), 21);
  StyleGraph(gRatioC.get(), 24);

  {
    TCanvas c("c_ratio_vs_run", "", 1600, 900);
    c.SetGrid();
    gRatioC->SetTitle("PHGarfield / CDB drift velocity;run;v_{PHGarfield}/v_{CDB}");
    gRatioC->Draw("APL");
    gRatioS->Draw("PL SAME");
    gRatioN->Draw("PL SAME");

    if (gRatioC->GetN() > 0)
    {
      const double xmin = gRatioC->GetX()[0];
      const double xmax = gRatioC->GetX()[gRatioC->GetN()-1];
      TLine one(xmin, 1.0, xmax, 1.0);
      one.SetLineStyle(2);
      one.Draw();
    }

    TLegend leg(0.14, 0.74, 0.35, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(gRatioC.get(), "combined", "lp");
    leg.AddEntry(gRatioS.get(), "south", "lp");
    leg.AddEntry(gRatioN.get(), "north", "lp");
    leg.Draw();

    SaveCanvas(c, plotDir, "02_garfield_over_cdb_vs_run");
  }

  auto gDiffS = MakeRelativeDifferenceRunGraph(runs, &RunData::garfieldSouth_cm_us, "g_diff_south");
  auto gDiffN = MakeRelativeDifferenceRunGraph(runs, &RunData::garfieldNorth_cm_us, "g_diff_north");
  auto gDiffC = MakeRelativeDifferenceRunGraph(runs, &RunData::garfieldCombined_cm_us, "g_diff_combined");
  StyleGraph(gDiffS.get(), 20);
  StyleGraph(gDiffN.get(), 21);
  StyleGraph(gDiffC.get(), 24);

  {
    TCanvas c("c_diff_vs_run", "", 1600, 900);
    c.SetGrid();
    gDiffC->SetTitle("PHGarfield - CDB relative difference;run;(v_{PHGarfield}-v_{CDB})/v_{CDB} [%]");
    gDiffC->Draw("APL");
    gDiffS->Draw("PL SAME");
    gDiffN->Draw("PL SAME");

    TLegend leg(0.14, 0.74, 0.35, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(gDiffC.get(), "combined", "lp");
    leg.AddEntry(gDiffS.get(), "south", "lp");
    leg.AddEntry(gDiffN.get(), "north", "lp");
    leg.Draw();

    SaveCanvas(c, plotDir, "03_relative_difference_vs_run");
  }

  // ----------------------------------------------------------------------
  // 3. Operating conditions / scaled field
  // ----------------------------------------------------------------------
  auto gField = MakeRunGraph(runs, &RunData::cmField_Vcm, "g_cmfield");
  auto gPressure = MakeRunGraph(runs, &RunData::pressure, "g_pressure");
  auto gTemp = MakeRunGraph(runs, &RunData::temperatureK, "g_temperatureK");
  StyleGraph(gField.get(), 20);
  StyleGraph(gPressure.get(), 20);
  StyleGraph(gTemp.get(), 20);

  {
    TCanvas c("c_field", "", 1600, 800);
    c.SetGrid();
    gField->SetTitle("Townsend-scaled nominal TPC field;run;E [V/cm]");
    gField->Draw("APL");
    SaveCanvas(c, plotDir, "04_scaled_field_vs_run");
  }
  {
    TCanvas c("c_pressure", "", 1600, 800);
    c.SetGrid();
    gPressure->SetTitle("TPC gas pressure from CDB;run;pressure [CDB units]");
    gPressure->Draw("APL");
    SaveCanvas(c, plotDir, "05_pressure_vs_run");
  }
  {
    TCanvas c("c_temperature", "", 1600, 800);
    c.SetGrid();
    gTemp->SetTitle("TPC gas temperature from CDB;run;T [K]");
    gTemp->Draw("APL");
    SaveCanvas(c, plotDir, "06_temperature_vs_run");
  }

  // ----------------------------------------------------------------------
  // 4. Average radial dependence over all runs, normalized to each run's CDB DV
  // ----------------------------------------------------------------------
  auto gR_S = AverageProfileRatio(runs, &RunData::vzRSouth, "g_avg_ratio_r_south");
  auto gR_N = AverageProfileRatio(runs, &RunData::vzRNorth, "g_avg_ratio_r_north");
  auto gTrajR_S = AverageProfileRatio(runs, &RunData::trajRSouth, "g_avg_traj_ratio_r_south");
  auto gTrajR_N = AverageProfileRatio(runs, &RunData::trajRNorth, "g_avg_traj_ratio_r_north");

  if (gR_S && gR_N)
  {
    StyleGraph(gR_S.get(), 20);
    StyleGraph(gR_N.get(), 21);
    TCanvas c("c_r_ratio", "", 1200, 850);
    c.SetGrid();
    gR_S->SetTitle("Mean local v_{z}/CDB vs radius across runs;r [cm];<v_{z,out}/v_{CDB}>");
    gR_S->Draw("AP");
    gR_N->Draw("P SAME");

    TLine one(18.0, 1.0, 82.0, 1.0);
    one.SetLineStyle(2);
    one.Draw();

    TLegend leg(0.15, 0.78, 0.32, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(gR_S.get(), "south", "p");
    leg.AddEntry(gR_N.get(), "north", "p");
    leg.Draw();

    SaveCanvas(c, plotDir, "07_average_local_ratio_vs_r");
  }

  if (gTrajR_S && gTrajR_N)
  {
    StyleGraph(gTrajR_S.get(), 20);
    StyleGraph(gTrajR_N.get(), 21);
    TCanvas c("c_traj_r_ratio", "", 1200, 850);
    c.SetGrid();
    gTrajR_S->SetTitle("CM-to-pad trajectory-average v_{z}/CDB vs pad radius;r_{pad} [cm];<v_{traj}/v_{CDB}>");
    gTrajR_S->Draw("AP");
    gTrajR_N->Draw("P SAME");

    TLine one(18.0, 1.0, 82.0, 1.0);
    one.SetLineStyle(2);
    one.Draw();

    TLegend leg(0.15, 0.78, 0.32, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(gTrajR_S.get(), "south", "p");
    leg.AddEntry(gTrajR_N.get(), "north", "p");
    leg.Draw();

    SaveCanvas(c, plotDir, "08_average_trajectory_ratio_vs_r");
  }

  // ----------------------------------------------------------------------
  // 5. Average z dependence over all runs
  // ----------------------------------------------------------------------
  auto gZ_S = AverageProfileRatio(runs, &RunData::vzZSouth, "g_avg_ratio_z_south");
  auto gZ_N = AverageProfileRatio(runs, &RunData::vzZNorth, "g_avg_ratio_z_north");

  if (gZ_S && gZ_N)
  {
    StyleGraph(gZ_S.get(), 20);
    StyleGraph(gZ_N.get(), 21);
    TCanvas c("c_z_ratio", "", 1200, 850);
    c.SetGrid();
    gZ_S->SetTitle("Mean local v_{z}/CDB vs z across runs;z [cm];<v_{z,out}/v_{CDB}>");
    gZ_S->Draw("AP");
    gZ_N->Draw("P SAME");

    TLine one(-120.0, 1.0, 120.0, 1.0);
    one.SetLineStyle(2);
    one.Draw();

    TLegend leg(0.15, 0.78, 0.32, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(gZ_S.get(), "south", "p");
    leg.AddEntry(gZ_N.get(), "north", "p");
    leg.Draw();

    SaveCanvas(c, plotDir, "09_average_local_ratio_vs_z");
  }

  // ----------------------------------------------------------------------
  // 6. Run-by-run maps: ratio vs radius and z
  // ----------------------------------------------------------------------
  auto hRunR_S = RunProfileRatioMap(
      runs, &RunData::vzRSouth,
      "h_run_r_ratio_south",
      "South: local v_{z}/CDB by run and radius;run;r [cm]",
      "r [cm]");
  auto hRunR_N = RunProfileRatioMap(
      runs, &RunData::vzRNorth,
      "h_run_r_ratio_north",
      "North: local v_{z}/CDB by run and radius;run;r [cm]",
      "r [cm]");
  auto hRunZ_S = RunProfileRatioMap(
      runs, &RunData::vzZSouth,
      "h_run_z_ratio_south",
      "South: local v_{z}/CDB by run and z;run;z [cm]",
      "z [cm]");
  auto hRunZ_N = RunProfileRatioMap(
      runs, &RunData::vzZNorth,
      "h_run_z_ratio_north",
      "North: local v_{z}/CDB by run and z;run;z [cm]",
      "z [cm]");

  auto drawMap = [&](TH2D *h, const char *canvasName, const std::string &base)
  {
    if (!h) return;
    const int width = std::max(1400, 20 * h->GetNbinsX());
    TCanvas c(canvasName, "", width, 900);
    c.SetRightMargin(0.14);
    c.SetBottomMargin(0.18);
    h->Draw("COLZ");
    SaveCanvas(c, plotDir, base);
  };

  drawMap(hRunR_S.get(), "c_run_r_s", "10_run_vs_r_ratio_south");
  drawMap(hRunR_N.get(), "c_run_r_n", "11_run_vs_r_ratio_north");
  drawMap(hRunZ_S.get(), "c_run_z_s", "12_run_vs_z_ratio_south");
  drawMap(hRunZ_N.get(), "c_run_z_n", "13_run_vs_z_ratio_north");

  // ----------------------------------------------------------------------
  // 7. Average 2D spatial ratio maps across all runs
  // ----------------------------------------------------------------------
  auto hRZ_S = AverageProfile2DRatio(
      runs, &RunData::vzRZSouth,
      "h_avg_rz_ratio_south",
      "South: mean local v_{z}/CDB;r [cm];z [cm]");
  auto hRZ_N = AverageProfile2DRatio(
      runs, &RunData::vzRZNorth,
      "h_avg_rz_ratio_north",
      "North: mean local v_{z}/CDB;r [cm];z [cm]");
  auto hRPhi_S = AverageProfile2DRatio(
      runs, &RunData::vzRPhiSouth,
      "h_avg_rphi_ratio_south",
      "South: mean local v_{z}/CDB;r [cm];#phi [rad]");
  auto hRPhi_N = AverageProfile2DRatio(
      runs, &RunData::vzRPhiNorth,
      "h_avg_rphi_ratio_north",
      "North: mean local v_{z}/CDB;r [cm];#phi [rad]");

  auto drawAverageMap = [&](TH2D *h, const char *canvasName, const std::string &base)
  {
    if (!h) return;
    TCanvas c(canvasName, "", 1200, 900);
    c.SetRightMargin(0.15);
    h->GetZaxis()->SetTitle("<v_{z,out}/v_{CDB}>");
    h->Draw("COLZ");
    SaveCanvas(c, plotDir, base);
  };

  drawAverageMap(hRZ_S.get(), "c_avg_rz_s", "14_average_rz_ratio_south");
  drawAverageMap(hRZ_N.get(), "c_avg_rz_n", "15_average_rz_ratio_north");
  drawAverageMap(hRPhi_S.get(), "c_avg_rphi_s", "16_average_rphi_ratio_south");
  drawAverageMap(hRPhi_N.get(), "c_avg_rphi_n", "17_average_rphi_ratio_north");

  // ----------------------------------------------------------------------
  // 8. Distribution of run-average PHGarfield/CDB ratios
  // ----------------------------------------------------------------------
  TH1D hRatioCombined("h_run_average_ratio_combined",
                      "Run-average PHGarfield/CDB ratio;v_{PHGarfield}/v_{CDB};runs",
                      160, 0.92, 1.08);
  TH1D hRatioSouth("h_run_average_ratio_south",
                   "Run-average PHGarfield/CDB ratio;v_{PHGarfield}/v_{CDB};runs",
                   160, 0.92, 1.08);
  TH1D hRatioNorth("h_run_average_ratio_north",
                   "Run-average PHGarfield/CDB ratio;v_{PHGarfield}/v_{CDB};runs",
                   160, 0.92, 1.08);

  for (const auto &r : runs)
  {
    const double cdb = 1000.0 * r.cdbDV_cm_ns;
    if (!finitePositive(cdb)) continue;
    if (std::isfinite(r.garfieldCombined_cm_us)) hRatioCombined.Fill(r.garfieldCombined_cm_us / cdb);
    if (std::isfinite(r.garfieldSouth_cm_us)) hRatioSouth.Fill(r.garfieldSouth_cm_us / cdb);
    if (std::isfinite(r.garfieldNorth_cm_us)) hRatioNorth.Fill(r.garfieldNorth_cm_us / cdb);
  }

  {
    TCanvas c("c_ratio_dist", "", 1100, 800);
    hRatioCombined.SetLineWidth(2);
    hRatioSouth.SetLineWidth(2);
    hRatioNorth.SetLineWidth(2);
    hRatioCombined.Draw("HIST");
    hRatioSouth.Draw("HIST SAME");
    hRatioNorth.Draw("HIST SAME");

    TLegend leg(0.68, 0.74, 0.88, 0.89);
    leg.SetBorderSize(0);
    leg.AddEntry(&hRatioCombined, "combined", "l");
    leg.AddEntry(&hRatioSouth, "south", "l");
    leg.AddEntry(&hRatioNorth, "north", "l");
    leg.Draw();

    SaveCanvas(c, plotDir, "18_run_average_ratio_distribution");
  }

  // Save all summary objects too.
  const std::string summaryRoot = plotDir + "/tpc_dv_postprocess.root";
  TFile fout(summaryRoot.c_str(), "RECREATE");

  gSouth->Write(); gNorth->Write(); gCombined->Write(); gCdb->Write(); gCdbT0->Write();
  gRatioS->Write(); gRatioN->Write(); gRatioC->Write();
  gDiffS->Write(); gDiffN->Write(); gDiffC->Write();
  gField->Write(); gPressure->Write(); gTemp->Write();

  if (gR_S) gR_S->Write();
  if (gR_N) gR_N->Write();
  if (gTrajR_S) gTrajR_S->Write();
  if (gTrajR_N) gTrajR_N->Write();
  if (gZ_S) gZ_S->Write();
  if (gZ_N) gZ_N->Write();

  if (hRunR_S) hRunR_S->Write();
  if (hRunR_N) hRunR_N->Write();
  if (hRunZ_S) hRunZ_S->Write();
  if (hRunZ_N) hRunZ_N->Write();
  if (hRZ_S) hRZ_S->Write();
  if (hRZ_N) hRZ_N->Write();
  if (hRPhi_S) hRPhi_S->Write();
  if (hRPhi_N) hRPhi_N->Write();

  hRatioCombined.Write();
  hRatioSouth.Write();
  hRatioNorth.Write();

  fout.Close();

  std::cout << "Plots saved in: " << plotDir << std::endl;
  std::cout << "Summary ROOT file: " << summaryRoot << std::endl;
}
