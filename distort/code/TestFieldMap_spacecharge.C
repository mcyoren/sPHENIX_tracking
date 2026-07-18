#include <QA.C>

#include <inttcalib/InttCalib.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

#include <GlobalVariables.C>
#include </sphenix/user/mitrankov/garf/PHGarfield/PHGarfield.h>

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TSystem.h>
#include <TVector3.h>

#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allutils.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffarawmodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(/sphenix/user/mitrankov/garf/PHGarfield/install/lib/libPHGarfield.so)

#define Nebdc 24
#define Nserver 2

TPolyLine3D* npoly3[48] = {};
TPolyLine3D* spoly3[48] = {};
TPolyLine* npoly2[48] = {};
TPolyLine* spoly2[48] = {};

TCanvas* canny = nullptr;
TCanvas* canny2 = nullptr;
TBox* boxer1 = nullptr;
TBox* boxer2 = nullptr;


// Find the existing stored point whose z is closest to target_z.
// No interpolation.
bool GetClosestPoint(const TPolyLine3D* line, double target_z,
                     double& r, double& dphi, double& actual_z,
                     int& selected_index, int& pad_index)
{
  if (!line || line->GetN() < 1 || !line->GetP()) return false;

  const int n = line->GetN();
  const float* p = line->GetP();

  const double z0 = p[2];
  const double z1 = p[3 * (n - 1) + 2];
  pad_index = std::abs(z0) >= std::abs(z1) ? 0 : n - 1;

  const double phi_pad = std::atan2(p[3 * pad_index + 1], p[3 * pad_index]);
  selected_index = -1;
  double best_dz = std::numeric_limits<double>::max();

  for (int i = 0; i < n; ++i)
  {
    const double dz = std::abs(p[3 * i + 2] - target_z);
    if (dz < best_dz)
    {
      best_dz = dz;
      selected_index = i;
    }
  }

  if (selected_index < 0) return false;

  const double x = p[3 * selected_index];
  const double y = p[3 * selected_index + 1];

  actual_z = p[3 * selected_index + 2];
  r = std::hypot(x, y);
  dphi = std::remainder(std::atan2(y, x) - phi_pad, 2.0 * TMath::Pi())*r;

  return true;
}


TPolyLine3D* MakeCircle3D(double r, double z, int color, int n = 180)
{
  auto* line = new TPolyLine3D(n + 1);
  line->SetLineColor(color);
  line->SetLineWidth(1);

  for (int i = 0; i <= n; ++i)
  {
    const double phi = 2.0 * TMath::Pi() * i / n;
    line->SetPoint(i, r * std::cos(phi), r * std::sin(phi), z);
  }

  return line;
}


TPolyLine3D* MakeZLine3D(double r, double phi, double zmin, double zmax, int color)
{
  auto* line = new TPolyLine3D(2);
  line->SetLineColor(color);
  line->SetLineWidth(1);
  line->SetPoint(0, r * std::cos(phi), r * std::sin(phi), zmin);
  line->SetPoint(1, r * std::cos(phi), r * std::sin(phi), zmax);
  return line;
}


void ConfigureMeasuredGraph(TGraph* graph, int color, int marker)
{
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(2.0);
  graph->SetLineWidth(2);
}


void AddDisplacementLine(std::vector<TLine*>& lines, double r_pad,
                         double r_measured, double dphi_cm, int color)
{
  auto* line = new TLine(r_pad, 0.0, r_measured, dphi_cm);
  line->SetLineColor(color);
  line->SetLineWidth(2);
  lines.push_back(line);
}


void UpdateCanvas(TCanvas* canvas)
{
  if (!canvas) return;
  canvas->cd();
  canvas->Modified();
  canvas->Update();
  gSystem->ProcessEvents();
}


void TestFieldMap_spacecharge(const double keff_side0 = 1.0,
                              const double keff_side1 = 1.0)
{
  recoConsts* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "FieldMapTest");
  rc->set_uint64Flag("TIMESTAMP", 1);

  auto* cdb = CDBInterface::instance();
  std::cout << "Field map URL:\n" << cdb->getUrl("FIELDMAP_TRACKING") << std::endl;

  Fun4AllServer* se = Fun4AllServer::instance();
  Enable::QA = false;
  Enable::CDB = true;

  // Input files.
  Fun4AllInputManager* in[Nebdc] = {};

  for (unsigned int ebdc = 0; ebdc < 2; ++ebdc)
  {
    for (unsigned int server = 0; server < 1; ++server)
    {
      const std::string input_name = std::format("ebdc{:02d}_{:1d}", ebdc, server);
      const std::string file_name = std::format(
          "DST_STREAMING_EVENT_ebdc{:02d}_{:1d}_"
          "run3auau_ana514_nocdbtag_v001-00075570-00000.root",
          ebdc, server);

      std::cout << file_name << " " << input_name << std::endl;

      in[ebdc] = new Fun4AllDstInputManager(input_name);
      in[ebdc]->fileopen(file_name);
      se->registerInputManager(in[ebdc]);
    }
  }

  // PHGarfield.
  const std::string electricFieldMap = "include/sphenix_rossegger_garfield_field.root";

  auto* phg = new PHGarfield("PHGarfield", electricFieldMap, keff_side0, keff_side1);

  TVector3 Northxyz(-0.001, -0.001, 1123.109);   // mm
  TVector3 Southxyz(-3.354, -0.673, -1137.382);  // mm
  TVector3 center = 0.5 * (Northxyz + Southxyz);
  center *= 0.1;  // mm -> cm

  phg->MoveTpc(center.X(), center.Y(), center.Z());
  phg->RotateTpc(0.0, 0.001485, 0.0);
  phg->RotateTpc(0.000298, 0.0, 0.0);
  phg->SetCMVoltageDefault(370.0);

  se->registerSubsystem(phg);
  se->run(4);

  constexpr int nLayers = 48;
  constexpr double zPadNorth = 102.0;
  constexpr double zPadSouth = -102.0;
  constexpr double zTarget0 = 0.0;
  constexpr double zTargetNorth50 = 50.0;
  constexpr double zTargetSouth50 = -50.0;
  constexpr double rInner = 20.0;
  constexpr double rOuter = 80.0;

  const std::string rootFileName = "reverse_drift_Polyfile_noTube.root";
  std::unique_ptr<TFile> output(TFile::Open(rootFileName.c_str(), "RECREATE"));

  if (!output || output->IsZombie())
  {
    std::cerr << "Could not create " << rootFileName << std::endl;
    return;
  }

  // ===========================================================================
  // Canvas 1: 3D
  // ===========================================================================

  canny = new TCanvas("canny", "Reverse-drift trajectories in 3D", 3000, 2500);
  canny->cd();
  canny->SetLeftMargin(0.08);
  canny->SetRightMargin(0.04);

  auto* frame3 = new TH3D("frame3",
      "Reverse-drift trajectories;x [cm];y [cm];z [cm]",
      10, -85.0, 85.0, 10, -85.0, 85.0, 10, -110.0, 110.0);

  frame3->SetDirectory(nullptr);
  frame3->SetStats(false);
  frame3->Draw();

  gPad->SetTheta(20.0);
  gPad->SetPhi(-45.0);

  std::vector<TPolyLine3D*> detectorLines = {
      MakeCircle3D(rInner, zPadSouth, kGray + 1),
      MakeCircle3D(rInner, zPadNorth, kGray + 1),
      MakeCircle3D(rOuter, zPadSouth, kGray + 1),
      MakeCircle3D(rOuter, zPadNorth, kGray + 1),
      MakeCircle3D(rInner, 0.0, kGray + 2),
      MakeCircle3D(rOuter, 0.0, kGray + 2)};

  for (int i = 0; i < 12; ++i)
  {
    const double phi = 2.0 * TMath::Pi() * i / 12.0;
    detectorLines.push_back(MakeZLine3D(rInner, phi, zPadSouth, zPadNorth, kGray + 1));
    detectorLines.push_back(MakeZLine3D(rOuter, phi, zPadSouth, zPadNorth, kGray + 1));
  }

  for (auto* line : detectorLines) line->Draw("same");

  // ===========================================================================
  // Canvas 2: r-z
  // ===========================================================================

  canny2 = new TCanvas("canny2", "Reverse-drift trajectories: r versus z", 3000, 2500);
  canny2->cd();
  canny2->SetGrid();

  auto* frameRZ = new TH2D("frameRZ",
      "Reverse-drift trajectories;z [cm];r [cm]",
      220, -110.0, 110.0, 150, 15.0, 85.0);

  frameRZ->SetDirectory(nullptr);
  frameRZ->SetStats(false);
  frameRZ->Draw();

  boxer1 = new TBox(-102.0, 20.0, 102.0, 80.0);
  boxer1->SetFillStyle(0);
  boxer1->SetLineColor(kGray + 1);
  boxer1->SetLineStyle(2);
  boxer1->Draw("same");

  boxer2 = new TBox(-0.2, 20.0, 0.2, 80.0);
  boxer2->SetFillStyle(0);
  boxer2->SetLineColor(kGray + 2);
  boxer2->SetLineStyle(2);
  boxer2->Draw("same");

  // ===========================================================================
  // Phi graphs
  // ===========================================================================

  auto* graphPhiZ0North = new TGraph();
  auto* graphPhiZ0South = new TGraph();
  auto* graphPhiZ50North = new TGraph();
  auto* graphPhiZ50South = new TGraph();

  graphPhiZ0North->SetName("graph_dphi_vs_r_z0_north");
  graphPhiZ0South->SetName("graph_dphi_vs_r_z0_south");
  graphPhiZ50North->SetName("graph_dphi_vs_r_zplus50_north");
  graphPhiZ50South->SetName("graph_dphi_vs_r_zminus50_south");

  ConfigureMeasuredGraph(graphPhiZ0North, kRed, 20);
  ConfigureMeasuredGraph(graphPhiZ0South, kCyan + 2, 21);
  ConfigureMeasuredGraph(graphPhiZ50North, kRed, 20);
  ConfigureMeasuredGraph(graphPhiZ50South, kCyan + 2, 21);

  // Pad-plane radii at Delta phi = 0.
  auto* graphPadRadiusZ0 = new TGraph();
  auto* graphPadRadiusZ50 = new TGraph();

  graphPadRadiusZ0->SetName("graph_pad_radius_z0");
  graphPadRadiusZ50->SetName("graph_pad_radius_z50");

  const int transparentGray = TColor::GetColorTransparent(kGray + 2, 0.60);
  const int transparenterGray = TColor::GetColorTransparent(kGray + 2, 0.20);

  for (auto* graph : {graphPadRadiusZ0, graphPadRadiusZ50})
  {
    graph->SetMarkerColor(transparentGray);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(2.5);
  }

  std::vector<TLine*> displacementLinesZ0;
  std::vector<TLine*> displacementLinesZ50;

  int nZ0North = 0, nZ0South = 0, nZ50North = 0, nZ50South = 0;

  // ===========================================================================
  // Generate drift lines
  // ===========================================================================

  for (int layer = 0; layer < nLayers; ++layer)
  {
    const double rPad = phg->GetRadius(layer);

    // One gray pad-radius marker per layer.
    graphPadRadiusZ0->SetPoint(layer, rPad, 0.0);
    graphPadRadiusZ50->SetPoint(layer, rPad, 0.0);

    // -------------------------------------------------------------------------
    // North
    // -------------------------------------------------------------------------

    npoly3[layer] = phg->ReverseDrift(0.0, rPad, zPadNorth);

    if (npoly3[layer] && npoly3[layer]->GetN() > 0 && npoly3[layer]->GetP())
    {
      npoly3[layer]->SetLineColor(kRed);
      npoly3[layer]->SetLineWidth(2);

      canny->cd();
      npoly3[layer]->Draw("same");

      const int n = npoly3[layer]->GetN();
      const float* p = npoly3[layer]->GetP();
      std::vector<double> z(n), r(n);

      for (int i = 0; i < n; ++i)
      {
        z[i] = p[3 * i + 2];
        r[i] = std::hypot(p[3 * i], p[3 * i + 1]);
      }

      npoly2[layer] = new TPolyLine(n, z.data(), r.data());
      npoly2[layer]->SetLineColor(kRed);
      npoly2[layer]->SetLineWidth(2);

      canny2->cd();
      npoly2[layer]->Draw("L same");

      double rMeasured = 0.0, dphi = 0.0, actualZ = 0.0;
      int selected = -1, pad = -1;

      if (GetClosestPoint(npoly3[layer], zTarget0, rMeasured, dphi, actualZ, selected, pad))
      {
        const double dphicm = 1. * dphi;
        graphPhiZ0North->SetPoint(nZ0North++, rMeasured, dphicm);
        AddDisplacementLine(displacementLinesZ0, rPad, rMeasured, dphicm, transparenterGray);
      }

      if (GetClosestPoint(npoly3[layer], zTargetNorth50, rMeasured, dphi, actualZ, selected, pad))
      {
        const double dphicm = 1. * dphi;
        graphPhiZ50North->SetPoint(nZ50North++, rMeasured, dphicm);
        AddDisplacementLine(displacementLinesZ50, rPad, rMeasured, dphicm, transparenterGray);
      }
    }

    // -------------------------------------------------------------------------
    // South
    // -------------------------------------------------------------------------

    spoly3[layer] = phg->ReverseDrift(0.0, rPad, zPadSouth);

    if (spoly3[layer] && spoly3[layer]->GetN() > 0 && spoly3[layer]->GetP())
    {
      spoly3[layer]->SetLineColor(kCyan);
      spoly3[layer]->SetLineWidth(2);

      canny->cd();
      spoly3[layer]->Draw("same");

      const int n = spoly3[layer]->GetN();
      const float* p = spoly3[layer]->GetP();
      std::vector<double> z(n), r(n);

      for (int i = 0; i < n; ++i)
      {
        z[i] = p[3 * i + 2];
        r[i] = std::hypot(p[3 * i], p[3 * i + 1]);
      }

      spoly2[layer] = new TPolyLine(n, z.data(), r.data());
      spoly2[layer]->SetLineColor(kCyan);
      spoly2[layer]->SetLineWidth(2);

      canny2->cd();
      spoly2[layer]->Draw("L same");

      double rMeasured = 0.0, dphi = 0.0, actualZ = 0.0;
      int selected = -1, pad = -1;

      if (GetClosestPoint(spoly3[layer], zTarget0, rMeasured, dphi, actualZ, selected, pad))
      {
        const double dphicm = 1.0 * dphi;
        graphPhiZ0South->SetPoint(nZ0South++, rMeasured, dphicm);
        AddDisplacementLine(displacementLinesZ0, rPad, rMeasured, dphicm, transparenterGray);
      }

      if (GetClosestPoint(spoly3[layer], zTargetSouth50, rMeasured, dphi, actualZ, selected, pad))
      {
        const double dphicm = 1.0 * dphi;
        graphPhiZ50South->SetPoint(nZ50South++, rMeasured, dphicm);
        AddDisplacementLine(displacementLinesZ50, rPad, rMeasured, dphicm, transparenterGray);
      }
    }
  }

  // ===========================================================================
  // Canvas 3: Delta phi near z = 0
  // ===========================================================================

  auto* canvasPhiZ0 = new TCanvas("canvasPhiZ0", "Phi distortion near z = 0", 2200, 1700);
  canvasPhiZ0->cd();
  canvasPhiZ0->SetGrid();

  auto* framePhiZ0 = new TH2D("framePhiZ0",
      "#Delta#phi and radial displacement near z = 0;"
      "r [cm];r#Delta#phi from pad plane [cm]",
      120, 20.0, 80.0, 400, -3.0, 3.0);

  framePhiZ0->SetDirectory(nullptr);
  framePhiZ0->SetStats(false);
  framePhiZ0->Draw();

  for (auto* line : displacementLinesZ0) line->Draw("same");

  if (graphPadRadiusZ0->GetN()) graphPadRadiusZ0->Draw("P same");
  if (graphPhiZ0North->GetN()) graphPhiZ0North->Draw("P same");
  if (graphPhiZ0South->GetN()) graphPhiZ0South->Draw("P same");

  auto* zeroLineZ0 = new TLine(20.0, 0.0, 80.0, 0.0);
  zeroLineZ0->SetLineColor(kGray + 2);
  zeroLineZ0->SetLineStyle(2);
  zeroLineZ0->SetLineWidth(2);
  zeroLineZ0->Draw("same");

  auto* legendPhiZ0 = new TLegend(0.52, 0.20, 0.89, 0.39);
  legendPhiZ0->SetBorderSize(0);
  legendPhiZ0->SetFillStyle(0);
  legendPhiZ0->AddEntry(graphPadRadiusZ0, "Pad-plane radius, #Delta#phi = 0", "p");
  legendPhiZ0->AddEntry(graphPhiZ0North, "North measured point", "p");
  legendPhiZ0->AddEntry(graphPhiZ0South, "South measured point", "p");
  legendPhiZ0->Draw();

  // ===========================================================================
  // Canvas 4: Delta phi near |z| = 50
  // ===========================================================================

  auto* canvasPhiZ50 = new TCanvas("canvasPhiZ50", "Phi distortion near |z| = 50", 2200, 1700);
  canvasPhiZ50->cd();
  canvasPhiZ50->SetGrid();

  auto* framePhiZ50 = new TH2D("framePhiZ50",
      "#Delta#phi and radial displacement near |z| = 50;"
      "r [cm];r#Delta#phi from pad plane [cm]",
      120, 20.0, 80.0, 400, -3.0, 3.0);

  framePhiZ50->SetDirectory(nullptr);
  framePhiZ50->SetStats(false);
  framePhiZ50->Draw();

  for (auto* line : displacementLinesZ50) line->Draw("same");

  if (graphPadRadiusZ50->GetN()) graphPadRadiusZ50->Draw("P same");
  if (graphPhiZ50North->GetN()) graphPhiZ50North->Draw("P same");
  if (graphPhiZ50South->GetN()) graphPhiZ50South->Draw("P same");

  auto* zeroLineZ50 = new TLine(20.0, 0.0, 80.0, 0.0);
  zeroLineZ50->SetLineColor(kGray + 2);
  zeroLineZ50->SetLineStyle(2);
  zeroLineZ50->SetLineWidth(2);
  zeroLineZ50->Draw("same");

  auto* legendPhiZ50 = new TLegend(0.52, 0.20, 0.89, 0.39);
  legendPhiZ50->SetBorderSize(0);
  legendPhiZ50->SetFillStyle(0);
  legendPhiZ50->AddEntry(graphPadRadiusZ50, "Pad-plane radius, #Delta#phi = 0", "p");
  legendPhiZ50->AddEntry(graphPhiZ50North, "North measured point", "p");
  legendPhiZ50->AddEntry(graphPhiZ50South, "South measured point", "p");
  legendPhiZ50->Draw();

  // ===========================================================================
  // Save
  // ===========================================================================

  UpdateCanvas(canny);
  UpdateCanvas(canny2);
  UpdateCanvas(canvasPhiZ0);
  UpdateCanvas(canvasPhiZ50);

  const std::string pdf3D = Form("PL_3D_V370_keff_%g_%g.pdf", keff_side0, keff_side1);
  const std::string pdfRZ = Form("PL_RZ_V370_keff_%g_%g.pdf", keff_side0, keff_side1);
  const std::string pdfPhiZ0 = Form("PL_dPhi_vs_r_z0_V370_keff_%g_%g.pdf", keff_side0, keff_side1);
  const std::string pdfPhiZ50 = Form("PL_dPhi_vs_r_z50_V370_keff_%g_%g.pdf", keff_side0, keff_side1);

  canny->SaveAs(pdf3D.c_str());
  canny2->SaveAs(pdfRZ.c_str());
  canvasPhiZ0->SaveAs(pdfPhiZ0.c_str());
  canvasPhiZ50->SaveAs(pdfPhiZ50.c_str());

  output->cd();

  canny->Write("canny", TObject::kOverwrite);
  canny2->Write("canny2", TObject::kOverwrite);
  canvasPhiZ0->Write("canvasPhiZ0", TObject::kOverwrite);
  canvasPhiZ50->Write("canvasPhiZ50", TObject::kOverwrite);

  graphPhiZ0North->Write("graph_dphi_vs_r_z0_north", TObject::kOverwrite);
  graphPhiZ0South->Write("graph_dphi_vs_r_z0_south", TObject::kOverwrite);
  graphPhiZ50North->Write("graph_dphi_vs_r_zplus50_north", TObject::kOverwrite);
  graphPhiZ50South->Write("graph_dphi_vs_r_zminus50_south", TObject::kOverwrite);
  graphPadRadiusZ0->Write("graph_pad_radius_z0", TObject::kOverwrite);
  graphPadRadiusZ50->Write("graph_pad_radius_z50", TObject::kOverwrite);

  output->Write();
  output->Close();

  std::cout
      << "\nSaved ROOT file: " << rootFileName
      << "\nSaved PDFs:\n  " << pdf3D
      << "\n  " << pdfRZ
      << "\n  " << pdfPhiZ0
      << "\n  " << pdfPhiZ50
      << "\nGraph points:"
      << "\n  North z~0: " << graphPhiZ0North->GetN()
      << "\n  South z~0: " << graphPhiZ0South->GetN()
      << "\n  North z~50: " << graphPhiZ50North->GetN()
      << "\n  South z~-50: " << graphPhiZ50South->GetN()
      << std::endl;
}