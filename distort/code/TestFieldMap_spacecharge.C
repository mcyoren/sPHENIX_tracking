#include <QA.C>

#include <inttcalib/InttCalib.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

#include </sphenix/user/mitrankov/garf/PHGarfield/PHGarfield.h>
#include <GlobalVariables.C>

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
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>

#include <array>
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

namespace
{
  constexpr int kNumEbdc = 24;
  constexpr int kNumLayers = 48;

  constexpr double kZPadNorth = 102.0;
  constexpr double kZPadSouth = -102.0;
  constexpr double kZTargetCenter = 0.0;
  constexpr double kZTargetNorth50 = 50.0;
  constexpr double kZTargetSouth50 = -50.0;

  constexpr double kRInner = 20.0;
  constexpr double kROuter = 80.0;

  constexpr int kNorthColor = kRed;
  constexpr int kSouthColor = kCyan + 2;

  const std::array<double, 3> kRdPhiTargetRadii = {30.0, 50.0, 70.0};

  std::array<TPolyLine3D *, kNumLayers> gNorthPoly3D = {};
  std::array<TPolyLine3D *, kNumLayers> gSouthPoly3D = {};
  std::array<TPolyLine *, kNumLayers> gNorthPolyRZ = {};
  std::array<TPolyLine *, kNumLayers> gSouthPolyRZ = {};

  TCanvas *gCanvas3D = nullptr;
  TCanvas *gCanvasRZ = nullptr;

  int GetPadPlanePointIndex(const TPolyLine3D *line)
  {
    if (!line || line->GetN() < 1 || !line->GetP())
    {
      return -1;
    }

    const int nPoints = line->GetN();
    const float *points = line->GetP();

    const double firstAbsZ = std::abs(points[2]);
    const double lastAbsZ = std::abs(points[3 * (nPoints - 1) + 2]);

    return firstAbsZ >= lastAbsZ ? 0 : nPoints - 1;
  }

  double GetPadPlanePhi(const TPolyLine3D *line)
  {
    const int padIndex = GetPadPlanePointIndex(line);
    if (padIndex < 0)
    {
      return 0.0;
    }

    const float *points = line->GetP();
    return std::atan2(points[3 * padIndex + 1], points[3 * padIndex]);
  }

  // Return the stored trajectory point whose z is closest to targetZ.
  // No interpolation is performed.
  bool GetClosestPoint(const TPolyLine3D *line, const double targetZ, double &radius, double &rdphi, double &actualZ, int &selectedIndex)
  {
    if (!line || line->GetN() < 1 || !line->GetP())
    {
      return false;
    }

    const int nPoints = line->GetN();
    const float *points = line->GetP();
    const double phiPad = GetPadPlanePhi(line);

    selectedIndex = -1;
    double smallestDeltaZ = std::numeric_limits<double>::max();

    for (int pointIndex = 0; pointIndex < nPoints; ++pointIndex)
    {
      const double z = points[3 * pointIndex + 2];
      const double deltaZ = std::abs(z - targetZ);

      if (deltaZ < smallestDeltaZ)
      {
        smallestDeltaZ = deltaZ;
        selectedIndex = pointIndex;
      }
    }

    if (selectedIndex < 0)
    {
      return false;
    }

    const double x = points[3 * selectedIndex];
    const double y = points[3 * selectedIndex + 1];

    actualZ = points[3 * selectedIndex + 2];
    radius = std::hypot(x, y);

    const double phi = std::atan2(y, x);
    rdphi = radius * std::remainder(phi - phiPad, 2.0 * TMath::Pi());

    return true;
  }

  TPolyLine3D *MakeCircle3D(const double radius, const double z, const int color, const int nSegments = 180)
  {
    auto *line = new TPolyLine3D(nSegments + 1);
    line->SetLineColor(color);
    line->SetLineWidth(1);

    for (int pointIndex = 0; pointIndex <= nSegments; ++pointIndex)
    {
      const double phi = 2.0 * TMath::Pi() * pointIndex / nSegments;
      line->SetPoint(pointIndex, radius * std::cos(phi), radius * std::sin(phi), z);
    }

    return line;
  }

  TPolyLine3D *MakeZLine3D(const double radius, const double phi, const double zMin, const double zMax, const int color)
  {
    auto *line = new TPolyLine3D(2);
    line->SetLineColor(color);
    line->SetLineWidth(1);

    line->SetPoint(0, radius * std::cos(phi), radius * std::sin(phi), zMin);
    line->SetPoint(1, radius * std::cos(phi), radius * std::sin(phi), zMax);

    return line;
  }

  TPolyLine *MakeRZPolyline(const TPolyLine3D *line, const int color)
  {
    if (!line || line->GetN() < 1 || !line->GetP())
    {
      return nullptr;
    }

    const int nPoints = line->GetN();
    const float *points = line->GetP();

    std::vector<double> z(nPoints);
    std::vector<double> radius(nPoints);

    for (int pointIndex = 0; pointIndex < nPoints; ++pointIndex)
    {
      const double x = points[3 * pointIndex];
      const double y = points[3 * pointIndex + 1];

      z[pointIndex] = points[3 * pointIndex + 2];
      radius[pointIndex] = std::hypot(x, y);
    }

    auto *rzLine = new TPolyLine(nPoints, z.data(), radius.data());
    rzLine->SetLineColor(color);
    rzLine->SetLineWidth(2);

    return rzLine;
  }

  TGraph *MakeRdPhiVsZGraph(const TPolyLine3D *line, const std::string &name, const int color)
  {
    auto *graph = new TGraph();
    graph->SetName(name.c_str());
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetLineWidth(3);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.7);

    if (!line || line->GetN() < 1 || !line->GetP())
    {
      return graph;
    }

    const int nPoints = line->GetN();
    const float *points = line->GetP();
    const double phiPad = GetPadPlanePhi(line);

    for (int pointIndex = 0; pointIndex < nPoints; ++pointIndex)
    {
      const double x = points[3 * pointIndex];
      const double y = points[3 * pointIndex + 1];
      const double z = points[3 * pointIndex + 2];

      const double radius = std::hypot(x, y);
      const double phi = std::atan2(y, x);
      const double rdphi = radius * std::remainder(phi - phiPad, 2.0 * TMath::Pi());

      graph->SetPoint(graph->GetN(), z, rdphi);
    }

    return graph;
  }

  void ConfigureMeasuredGraph(TGraph *graph, const int color, const int markerStyle)
  {
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(2.0);
    graph->SetLineWidth(2);
  }

  void AddDisplacementLine(std::vector<TLine *> &lines, const double padRadius, const double measuredRadius, const double rdphi, const int color)
  {
    auto *line = new TLine(padRadius, 0.0, measuredRadius, rdphi);
    line->SetLineColor(color);
    line->SetLineWidth(2);
    lines.push_back(line);
  }

  void UpdateCanvas(TCanvas *canvas)
  {
    if (!canvas)
    {
      return;
    }

    canvas->cd();
    canvas->Modified();
    canvas->Update();
    gSystem->ProcessEvents();
  }

  int FindClosestLayer(PHGarfield *phg, const double targetRadius)
  {
    int closestLayer = -1;
    double smallestDifference = std::numeric_limits<double>::max();

    for (int layer = 0; layer < kNumLayers; ++layer)
    {
      const double radius = phg->GetRadius(layer);
      const double difference = std::abs(radius - targetRadius);

      if (difference < smallestDifference)
      {
        smallestDifference = difference;
        closestLayer = layer;
      }
    }

    return closestLayer;
  }

  void DrawTrajectory(TPolyLine3D *trajectory, TPolyLine *&rzTrajectory, TCanvas *canvas3D, TCanvas *canvasRZ, const int color)
  {
    if (!trajectory || trajectory->GetN() < 1 || !trajectory->GetP())
    {
      return;
    }

    trajectory->SetLineColor(color);
    trajectory->SetLineWidth(2);

    canvas3D->cd();
    trajectory->Draw("same");

    rzTrajectory = MakeRZPolyline(trajectory, color);
    if (rzTrajectory)
    {
      canvasRZ->cd();
      rzTrajectory->Draw("L same");
    }
  }

  void FillRdPhiAtZ(TPolyLine3D *trajectory, const double targetZ, const double padRadius, TGraph *graph, int &graphPoint, std::vector<TLine *> &displacementLines, const int displacementColor)
  {
    double measuredRadius = 0.0;
    double rdphi = 0.0;
    double actualZ = 0.0;
    int selectedIndex = -1;

    if (!GetClosestPoint(trajectory, targetZ, measuredRadius, rdphi, actualZ, selectedIndex))
    {
      return;
    }

    graph->SetPoint(graphPoint++, measuredRadius, rdphi);

    AddDisplacementLine(displacementLines, padRadius, measuredRadius, rdphi, displacementColor);
  }
} // namespace

void TestFieldMap_spacecharge(const double keff_side0 = 1.0, const double keff_side1 = 1.0)
{
  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "FieldMapTest");
  rc->set_uint64Flag("TIMESTAMP", 1);

  auto *cdb = CDBInterface::instance();
  std::cout << "Field map URL:\n"
            << cdb->getUrl("FIELDMAP_TRACKING") << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();

  Enable::QA = false;
  Enable::CDB = true;

  // ===========================================================================
  // Input
  // ===========================================================================

  std::array<Fun4AllInputManager *, kNumEbdc> inputManagers = {};

  for (unsigned int ebdc = 0; ebdc < 2; ++ebdc)
  {
    for (unsigned int server = 0; server < 1; ++server)
    {
      const std::string inputName = std::format("ebdc{:02d}_{:1d}", ebdc, server);

      const std::string fileName = std::format("DST_STREAMING_EVENT_ebdc{:02d}_{:1d}_"
                                               "run3auau_ana514_nocdbtag_v001-00075570-00000.root",
                                               ebdc, server);

      std::cout << fileName << " " << inputName << std::endl;

      inputManagers[ebdc] = new Fun4AllDstInputManager(inputName);
      inputManagers[ebdc]->fileopen(fileName);
      se->registerInputManager(inputManagers[ebdc]);
    }
  }

  // ===========================================================================
  // PHGarfield configuration
  // ===========================================================================

  const std::string electricFieldMap = "";
  //"include/sphenix_rossegger_garfield_field_2p0.root";
  //const std::string field3DSide0 = "include/sphenix_3d_ibf_field_side0_South_v3.root";
  //const std::string field3DSide1 = "include/sphenix_3d_ibf_field_side1_North_v3.root";
  const std::string field3DSide0 = "include/ibf_side0_3d.root";
  const std::string field3DSide1 = "include/ibf_side1_3d.root";

  auto *phg = new PHGarfield("PHGarfield", electricFieldMap, keff_side0, keff_side1);

  phg->SetElectricFieldMap3D(field3DSide0, field3DSide1);

  const std::string framefield3DSide0 = "include/frames_side0_3d_v1.root";
  const std::string framefield3DSide1 = "include/frames_side1_3d_v1.root";

  phg->SetFrameElectricFieldMap3D(framefield3DSide0, framefield3DSide1);

  phg->SetFrameChargeScale(-3000);

  TVector3 northPositionMm(-0.001, -0.001, 1123.109);
  TVector3 southPositionMm(-3.354, -0.673, -1137.382);

  TVector3 tpcCenterCm = 0.5 * (northPositionMm + southPositionMm);
  tpcCenterCm *= 0.1; // mm -> cm

  // phg->MoveTpc(tpcCenterCm.X(),tpcCenterCm.Y(),tpcCenterCm.Z());
  // phg->RotateTpc(0.0, 0.001485, 0.0);
  // phg->RotateTpc(0.000298, 0.0, 0.0);
  phg->SetCMVoltageDefault(380.0);

  se->registerSubsystem(phg);
  se->run(4);

  // ===========================================================================
  // Output file
  // ===========================================================================

  const std::string rootFileName = Form("reverse_drift_Polyfile_V370_keff_%g_%g.root", keff_side0, keff_side1);

  std::unique_ptr<TFile> output(TFile::Open(rootFileName.c_str(), "RECREATE"));

  if (!output || output->IsZombie())
  {
    std::cerr << "Could not create " << rootFileName << std::endl;
    return;
  }

  // ===========================================================================
  // Canvas 1: 3D trajectories
  // ===========================================================================

  gCanvas3D = new TCanvas("canvas3D", "Reverse-drift trajectories in 3D", 3000, 2500);

  gCanvas3D->cd();
  gCanvas3D->SetLeftMargin(0.08);
  gCanvas3D->SetRightMargin(0.04);

  auto *frame3D = new TH3D("frame3D", "Reverse-drift trajectories;x [cm];y [cm];z [cm]", 10, -85.0, 85.0, 10, -85.0, 85.0, 10, -110.0, 110.0);

  frame3D->SetDirectory(nullptr);
  frame3D->SetStats(false);
  frame3D->Draw();

  gPad->SetTheta(20.0);
  gPad->SetPhi(-45.0);

  std::vector<TPolyLine3D *> detectorLines = {MakeCircle3D(kRInner, kZPadSouth, kGray + 1), MakeCircle3D(kRInner, kZPadNorth, kGray + 1), MakeCircle3D(kROuter, kZPadSouth, kGray + 1),
                                              MakeCircle3D(kROuter, kZPadNorth, kGray + 1), MakeCircle3D(kRInner, 0.0, kGray + 2), MakeCircle3D(kROuter, 0.0, kGray + 2)};

  for (int sector = 0; sector < 12; ++sector)
  {
    const double phi = 2.0 * TMath::Pi() * sector / 12.0;

    detectorLines.push_back(MakeZLine3D(kRInner, phi, kZPadSouth, kZPadNorth, kGray + 1));

    detectorLines.push_back(MakeZLine3D(kROuter, phi, kZPadSouth, kZPadNorth, kGray + 1));
  }

  for (auto *line : detectorLines)
  {
    line->Draw("same");
  }

  // ===========================================================================
  // Canvas 2: r versus z
  // ===========================================================================

  gCanvasRZ = new TCanvas("canvasRZ", "Reverse-drift trajectories: r versus z", 3000, 2500);

  gCanvasRZ->cd();
  gCanvasRZ->SetGrid();

  auto *frameRZ = new TH2D("frameRZ", "Reverse-drift trajectories;z [cm];r [cm]", 220, -110.0, 110.0, 150, 15.0, 85.0);

  frameRZ->SetDirectory(nullptr);
  frameRZ->SetStats(false);
  frameRZ->Draw();

  auto *activeVolumeBox = new TBox(kZPadSouth, kRInner, kZPadNorth, kROuter);

  activeVolumeBox->SetFillStyle(0);
  activeVolumeBox->SetLineColor(kGray + 1);
  activeVolumeBox->SetLineStyle(2);
  activeVolumeBox->Draw("same");

  auto *centralMembraneBox = new TBox(-0.2, kRInner, 0.2, kROuter);

  centralMembraneBox->SetFillStyle(0);
  centralMembraneBox->SetLineColor(kGray + 2);
  centralMembraneBox->SetLineStyle(2);
  centralMembraneBox->Draw("same");

  // ===========================================================================
  // rDeltaPhi versus r at selected z positions
  // ===========================================================================

  auto *graphRdPhiZ0North = new TGraph();
  auto *graphRdPhiZ0South = new TGraph();
  auto *graphRdPhiZ50North = new TGraph();
  auto *graphRdPhiZ50South = new TGraph();

  graphRdPhiZ0North->SetName("graph_rdphi_vs_r_z0_north");
  graphRdPhiZ0South->SetName("graph_rdphi_vs_r_z0_south");
  graphRdPhiZ50North->SetName("graph_rdphi_vs_r_zplus50_north");
  graphRdPhiZ50South->SetName("graph_rdphi_vs_r_zminus50_south");

  ConfigureMeasuredGraph(graphRdPhiZ0North, kNorthColor, 20);
  ConfigureMeasuredGraph(graphRdPhiZ0South, kSouthColor, 21);
  ConfigureMeasuredGraph(graphRdPhiZ50North, kNorthColor, 20);
  ConfigureMeasuredGraph(graphRdPhiZ50South, kSouthColor, 21);

  auto *graphPadRadiusZ0 = new TGraph();
  auto *graphPadRadiusZ50 = new TGraph();

  graphPadRadiusZ0->SetName("graph_pad_radius_z0");
  graphPadRadiusZ50->SetName("graph_pad_radius_z50");

  const int transparentGray = TColor::GetColorTransparent(kGray + 2, 0.60);

  const int veryTransparentGray = TColor::GetColorTransparent(kGray + 2, 0.20);

  for (auto *graph : {graphPadRadiusZ0, graphPadRadiusZ50})
  {
    graph->SetMarkerColor(transparentGray);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(2.5);
  }

  std::vector<TLine *> displacementLinesZ0;
  std::vector<TLine *> displacementLinesZ50;

  int nZ0North = 0;
  int nZ0South = 0;
  int nZ50North = 0;
  int nZ50South = 0;

  // ===========================================================================
  // Generate all drift trajectories
  // ===========================================================================

  for (int layer = 0; layer < kNumLayers; ++layer)
  {
    const double padRadius = phg->GetRadius(layer);

    graphPadRadiusZ0->SetPoint(layer, padRadius, 0.0);
    graphPadRadiusZ50->SetPoint(layer, padRadius, 0.0);

    // North trajectory.
    gNorthPoly3D[layer] = phg->ReverseDrift(0.0, padRadius, kZPadNorth);

    DrawTrajectory(gNorthPoly3D[layer], gNorthPolyRZ[layer], gCanvas3D, gCanvasRZ, kNorthColor);

    FillRdPhiAtZ(gNorthPoly3D[layer], kZTargetCenter, padRadius, graphRdPhiZ0North, nZ0North, displacementLinesZ0, veryTransparentGray);

    FillRdPhiAtZ(gNorthPoly3D[layer], kZTargetNorth50, padRadius, graphRdPhiZ50North, nZ50North, displacementLinesZ50, veryTransparentGray);

    // South trajectory.
    gSouthPoly3D[layer] = phg->ReverseDrift(0.0, padRadius, kZPadSouth);

    DrawTrajectory(gSouthPoly3D[layer], gSouthPolyRZ[layer], gCanvas3D, gCanvasRZ, kSouthColor);

    FillRdPhiAtZ(gSouthPoly3D[layer], kZTargetCenter, padRadius, graphRdPhiZ0South, nZ0South, displacementLinesZ0, veryTransparentGray);

    FillRdPhiAtZ(gSouthPoly3D[layer], kZTargetSouth50, padRadius, graphRdPhiZ50South, nZ50South, displacementLinesZ50, veryTransparentGray);
  }

  // ===========================================================================
  // Canvas 3: rDeltaPhi versus r near z = 0
  // ===========================================================================

  auto *canvasRdPhiZ0 = new TCanvas("canvasRdPhiZ0", "rDeltaPhi distortion near z = 0", 2200, 1700);

  canvasRdPhiZ0->cd();
  canvasRdPhiZ0->SetGrid();

  auto *frameRdPhiZ0 = new TH2D("frameRdPhiZ0",
                                "Azimuthal and radial displacement near z = 0;"
                                "r [cm];r#Delta#phi from pad plane [cm]",
                                120, 20.0, 80.0, 400, -3.0, 3.0);

  frameRdPhiZ0->SetDirectory(nullptr);
  frameRdPhiZ0->SetStats(false);
  frameRdPhiZ0->Draw();

  for (auto *line : displacementLinesZ0)
  {
    line->Draw("same");
  }

  if (graphPadRadiusZ0->GetN() > 0)
  {
    graphPadRadiusZ0->Draw("P same");
  }

  if (graphRdPhiZ0North->GetN() > 0)
  {
    graphRdPhiZ0North->Draw("P same");
  }

  if (graphRdPhiZ0South->GetN() > 0)
  {
    graphRdPhiZ0South->Draw("P same");
  }

  auto *zeroLineZ0 = new TLine(20.0, 0.0, 80.0, 0.0);
  zeroLineZ0->SetLineColor(kGray + 2);
  zeroLineZ0->SetLineStyle(2);
  zeroLineZ0->SetLineWidth(2);
  zeroLineZ0->Draw("same");

  auto *legendRdPhiZ0 = new TLegend(0.52, 0.20, 0.89, 0.39);
  legendRdPhiZ0->SetBorderSize(0);
  legendRdPhiZ0->SetFillStyle(0);
  legendRdPhiZ0->AddEntry(graphPadRadiusZ0, "Pad-plane radius, r#Delta#phi = 0", "p");
  legendRdPhiZ0->AddEntry(graphRdPhiZ0North, "North measured point", "p");
  legendRdPhiZ0->AddEntry(graphRdPhiZ0South, "South measured point", "p");
  legendRdPhiZ0->Draw();

  // ===========================================================================
  // Canvas 4: rDeltaPhi versus r near |z| = 50 cm
  // ===========================================================================

  auto *canvasRdPhiZ50 = new TCanvas("canvasRdPhiZ50", "rDeltaPhi distortion near |z| = 50 cm", 2200, 1700);

  canvasRdPhiZ50->cd();
  canvasRdPhiZ50->SetGrid();

  auto *frameRdPhiZ50 = new TH2D("frameRdPhiZ50",
                                 "Azimuthal and radial displacement near |z| = 50 cm;"
                                 "r [cm];r#Delta#phi from pad plane [cm]",
                                 120, 20.0, 80.0, 400, -3.0, 3.0);

  frameRdPhiZ50->SetDirectory(nullptr);
  frameRdPhiZ50->SetStats(false);
  frameRdPhiZ50->Draw();

  for (auto *line : displacementLinesZ50)
  {
    line->Draw("same");
  }

  if (graphPadRadiusZ50->GetN() > 0)
  {
    graphPadRadiusZ50->Draw("P same");
  }

  if (graphRdPhiZ50North->GetN() > 0)
  {
    graphRdPhiZ50North->Draw("P same");
  }

  if (graphRdPhiZ50South->GetN() > 0)
  {
    graphRdPhiZ50South->Draw("P same");
  }

  auto *zeroLineZ50 = new TLine(20.0, 0.0, 80.0, 0.0);
  zeroLineZ50->SetLineColor(kGray + 2);
  zeroLineZ50->SetLineStyle(2);
  zeroLineZ50->SetLineWidth(2);
  zeroLineZ50->Draw("same");

  auto *legendRdPhiZ50 = new TLegend(0.52, 0.20, 0.89, 0.39);
  legendRdPhiZ50->SetBorderSize(0);
  legendRdPhiZ50->SetFillStyle(0);
  legendRdPhiZ50->AddEntry(graphPadRadiusZ50, "Pad-plane radius, r#Delta#phi = 0", "p");
  legendRdPhiZ50->AddEntry(graphRdPhiZ50North, "North measured point", "p");
  legendRdPhiZ50->AddEntry(graphRdPhiZ50South, "South measured point", "p");
  legendRdPhiZ50->Draw();

  // ===========================================================================
  // Canvas 5: rDeltaPhi versus z near r = 30, 50, and 70 cm
  // ===========================================================================

  std::array<int, 3> selectedLayers = {};
  std::array<double, 3> selectedRadii = {};
  std::array<TGraph *, 3> graphRdPhiVsZNorth = {};
  std::array<TGraph *, 3> graphRdPhiVsZSouth = {};

  for (std::size_t radiusIndex = 0; radiusIndex < kRdPhiTargetRadii.size(); ++radiusIndex)
  {
    selectedLayers[radiusIndex] = FindClosestLayer(phg, kRdPhiTargetRadii[radiusIndex]);

    selectedRadii[radiusIndex] = phg->GetRadius(selectedLayers[radiusIndex]);

    graphRdPhiVsZNorth[radiusIndex] = MakeRdPhiVsZGraph(gNorthPoly3D[selectedLayers[radiusIndex]], Form("graph_rdphi_vs_z_r%.1f_north", selectedRadii[radiusIndex]), kNorthColor);

    graphRdPhiVsZSouth[radiusIndex] = MakeRdPhiVsZGraph(gSouthPoly3D[selectedLayers[radiusIndex]], Form("graph_rdphi_vs_z_r%.1f_south", selectedRadii[radiusIndex]), kSouthColor);
  }

  auto *canvasRdPhiVsZ = new TCanvas("canvasRdPhiVsZ", "rDeltaPhi versus z at selected radii", 3000, 1000);

  canvasRdPhiVsZ->Divide(3, 1);

  std::array<TH2D *, 3> frameRdPhiVsZ = {};
  std::array<TLine *, 3> zeroLinesRdPhiVsZ = {};
  std::array<TLegend *, 3> legendsRdPhiVsZ = {};

  for (std::size_t radiusIndex = 0; radiusIndex < kRdPhiTargetRadii.size(); ++radiusIndex)
  {
    canvasRdPhiVsZ->cd(radiusIndex + 1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.10);

    frameRdPhiVsZ[radiusIndex] = new TH2D(Form("frameRdPhiVsZ_%zu", radiusIndex),
                                          Form("Layer %d, r_{pad} = %.3f cm;"
                                               "z [cm];r#Delta#phi from pad plane [cm]",
                                               selectedLayers[radiusIndex], selectedRadii[radiusIndex]),
                                          220, -110.0, 110.0, 400, -3.0, 3.0);

    frameRdPhiVsZ[radiusIndex]->SetDirectory(nullptr);
    frameRdPhiVsZ[radiusIndex]->SetStats(false);
    frameRdPhiVsZ[radiusIndex]->Draw();

    zeroLinesRdPhiVsZ[radiusIndex] = new TLine(-110.0, 0.0, 110.0, 0.0);

    zeroLinesRdPhiVsZ[radiusIndex]->SetLineColor(kGray + 2);
    zeroLinesRdPhiVsZ[radiusIndex]->SetLineStyle(2);
    zeroLinesRdPhiVsZ[radiusIndex]->SetLineWidth(2);
    zeroLinesRdPhiVsZ[radiusIndex]->Draw("same");

    if (graphRdPhiVsZNorth[radiusIndex]->GetN() > 0)
    {
      graphRdPhiVsZNorth[radiusIndex]->Draw("LP same");
    }

    if (graphRdPhiVsZSouth[radiusIndex]->GetN() > 0)
    {
      graphRdPhiVsZSouth[radiusIndex]->Draw("LP same");
    }

    legendsRdPhiVsZ[radiusIndex] = new TLegend(0.17, 0.73, 0.48, 0.89);

    legendsRdPhiVsZ[radiusIndex]->SetBorderSize(0);
    legendsRdPhiVsZ[radiusIndex]->SetFillStyle(0);
    legendsRdPhiVsZ[radiusIndex]->AddEntry(graphRdPhiVsZNorth[radiusIndex], "North", "lp");
    legendsRdPhiVsZ[radiusIndex]->AddEntry(graphRdPhiVsZSouth[radiusIndex], "South", "lp");
    legendsRdPhiVsZ[radiusIndex]->Draw();
  }

  // ===========================================================================
  // Save canvases and graphs
  // ===========================================================================

  UpdateCanvas(gCanvas3D);
  UpdateCanvas(gCanvasRZ);
  UpdateCanvas(canvasRdPhiZ0);
  UpdateCanvas(canvasRdPhiZ50);
  UpdateCanvas(canvasRdPhiVsZ);

  const std::string pdf3D = Form("PL_3D_V370_keff_%g_%g.pdf", keff_side0, keff_side1);

  const std::string pdfRZ = Form("PL_RZ_V370_keff_%g_%g.pdf", keff_side0, keff_side1);

  const std::string pdfRdPhiZ0 = Form("PL_rDPhi_vs_r_z0_V370_keff_%g_%g.pdf", keff_side0, keff_side1);

  const std::string pdfRdPhiZ50 = Form("PL_rDPhi_vs_r_z50_V370_keff_%g_%g.pdf", keff_side0, keff_side1);

  const std::string pdfRdPhiVsZ = Form("PL_rDPhi_vs_z_r30_r50_r70_V370_keff_%g_%g.pdf", keff_side0, keff_side1);

  gCanvas3D->SaveAs(pdf3D.c_str());
  gCanvasRZ->SaveAs(pdfRZ.c_str());
  canvasRdPhiZ0->SaveAs(pdfRdPhiZ0.c_str());
  canvasRdPhiZ50->SaveAs(pdfRdPhiZ50.c_str());
  canvasRdPhiVsZ->SaveAs(pdfRdPhiVsZ.c_str());

  output->cd();

  gCanvas3D->Write("canvas3D", TObject::kOverwrite);
  gCanvasRZ->Write("canvasRZ", TObject::kOverwrite);
  canvasRdPhiZ0->Write("canvasRdPhiZ0", TObject::kOverwrite);
  canvasRdPhiZ50->Write("canvasRdPhiZ50", TObject::kOverwrite);
  canvasRdPhiVsZ->Write("canvasRdPhiVsZ", TObject::kOverwrite);

  graphRdPhiZ0North->Write("graph_rdphi_vs_r_z0_north", TObject::kOverwrite);

  graphRdPhiZ0South->Write("graph_rdphi_vs_r_z0_south", TObject::kOverwrite);

  graphRdPhiZ50North->Write("graph_rdphi_vs_r_zplus50_north", TObject::kOverwrite);

  graphRdPhiZ50South->Write("graph_rdphi_vs_r_zminus50_south", TObject::kOverwrite);

  graphPadRadiusZ0->Write("graph_pad_radius_z0", TObject::kOverwrite);

  graphPadRadiusZ50->Write("graph_pad_radius_z50", TObject::kOverwrite);

  for (std::size_t radiusIndex = 0; radiusIndex < kRdPhiTargetRadii.size(); ++radiusIndex)
  {
    graphRdPhiVsZNorth[radiusIndex]->Write(graphRdPhiVsZNorth[radiusIndex]->GetName(), TObject::kOverwrite);

    graphRdPhiVsZSouth[radiusIndex]->Write(graphRdPhiVsZSouth[radiusIndex]->GetName(), TObject::kOverwrite);
  }

  output->Write();
  output->Close();

  std::cout << "\nSaved ROOT file: " << rootFileName << "\nSaved PDFs:"
            << "\n  " << pdf3D << "\n  " << pdfRZ << "\n  " << pdfRdPhiZ0 << "\n  " << pdfRdPhiZ50 << "\n  " << pdfRdPhiVsZ << "\n\nGraph points:"
            << "\n  North z ~ 0: " << graphRdPhiZ0North->GetN() << "\n  South z ~ 0: " << graphRdPhiZ0South->GetN() << "\n  North z ~ +50 cm: " << graphRdPhiZ50North->GetN()
            << "\n  South z ~ -50 cm: " << graphRdPhiZ50South->GetN() << "\n\nSelected layers for rDeltaPhi versus z:";

  for (std::size_t radiusIndex = 0; radiusIndex < kRdPhiTargetRadii.size(); ++radiusIndex)
  {
    std::cout << "\n  requested r ~ " << kRdPhiTargetRadii[radiusIndex] << " cm: layer " << selectedLayers[radiusIndex] << ", actual pad radius = " << selectedRadii[radiusIndex] << " cm";
  }

  std::cout << std::endl;
}
