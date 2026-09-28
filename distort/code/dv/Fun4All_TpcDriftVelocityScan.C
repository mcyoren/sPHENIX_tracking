// Fun4All_TpcDriftVelocityScan.C
//
// One-run PHGarfield drift-velocity scan.
//
// Goals:
//   * use PHGarfield's run-dependent CDB configuration
//   * set the nominal CM field with Townsend E/N scaling relative to run 79513
//   * compare PHGarfield drift velocity with TPC_DRIFT_VELOCITY from CDB
//   * read TPC_TZERO_OFFSET directly from CDB
//   * optionally correct the CDB drift velocity for a separately measured t0
//   * calculate local velocity from every ReverseDrift polyline step
//   * save 1D velocity distributions and <vz> maps vs r-z and r-phi
//
// Recommended usage: ONE ROOT/Fun4All process per run. See run_tpc_dv_scan.sh.

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <phfield/PHFieldConfig.h>
#include <phfield/PHFieldConfigv1.h>
#include <phfield/PHFieldUtility.h>
#include <phgarfield/PHGarfield.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TParameter.h>
#include <TPolyLine3D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libphfield.so)
R__LOAD_LIBRARY(libPHGarfield.so)

namespace TpcDVScan
{
  constexpr int kReferenceRun = 79513;
  constexpr double kReferenceCMField_Vcm = 377.0;
  constexpr int kNLayers = 48;

  struct RunConditions
  {
    bool valid{false};
    double pressure{std::numeric_limits<double>::quiet_NaN()};
    double temperatureRaw{std::numeric_limits<double>::quiet_NaN()};
    double temperatureK{std::numeric_limits<double>::quiet_NaN()};
    int nSamples{0};
  };

  struct CdbRecoConstants
  {
    double driftVelocity_cm_ns{std::numeric_limits<double>::quiet_NaN()};
    double tzero_ns{std::numeric_limits<double>::quiet_NaN()};
  };

  double ToKelvin(const double temperature)
  {
    // TPC_CONDITIONS currently does not document the temperature unit in its
    // public README. This makes the macro robust to the common cases:
    // values below 200 are interpreted as deg C, otherwise as K.
    if (!std::isfinite(temperature)) return temperature;
    return temperature < 200.0 ? temperature + 273.15 : temperature;
  }

  void SetCdbRun(const int runnumber, const std::string &globalTag)
  {
    auto *rc = recoConsts::instance();
    rc->set_StringFlag("CDB_GLOBALTAG", globalTag);
    rc->set_IntFlag("RUNNUMBER", runnumber);
    rc->set_uint64Flag("TIMESTAMP", static_cast<uint64_t>(runnumber));
  }

  RunConditions ReadAverageRunConditions(const int runnumber,
                                         const std::string &globalTag)
  {
    SetCdbRun(runnumber, globalTag);

    RunConditions result;
    const std::string url = CDBInterface::instance()->getUrl("TPC_CONDITIONS");
    if (url.empty())
    {
      std::cout << "TPC_CONDITIONS is not available for run " << runnumber << std::endl;
      return result;
    }

    CDBTTree tree(url);
    tree.LoadCalibrations();

    const std::size_t nEntries = tree.GetUInt64EntryMap().size();
    double sumP = 0.0;
    double sumT = 0.0;
    int n = 0;

    for (std::size_t i = 0; i < nEntries; ++i)
    {
      const double p = tree.GetFloatValue(static_cast<int>(i), "gas_pressure");
      const double t = tree.GetFloatValue(static_cast<int>(i), "gas_temperature");
      if (!std::isfinite(p) || !std::isfinite(t) || p <= 0.0)
      {
        continue;
      }
      sumP += p;
      sumT += t;
      ++n;
    }

    if (n == 0)
    {
      std::cout << "No valid P/T samples in TPC_CONDITIONS for run " << runnumber << std::endl;
      return result;
    }

    result.valid = true;
    result.pressure = sumP / n;
    result.temperatureRaw = sumT / n;
    result.temperatureK = ToKelvin(result.temperatureRaw);
    result.nSamples = n;
    return result;
  }

  CdbRecoConstants ReadCdbRecoConstants(const int runnumber,
                                        const std::string &globalTag)
  {
    SetCdbRun(runnumber, globalTag);

    CdbRecoConstants result;
    auto *cdb = CDBInterface::instance();

    const std::string dvUrl = cdb->getUrl("TPC_DRIFT_VELOCITY");
    if (!dvUrl.empty())
    {
      CDBTTree tree(dvUrl);
      tree.LoadCalibrations();
      result.driftVelocity_cm_ns = tree.GetSingleFloatValue("tpc_drift_velocity");
    }

    const std::string t0Url = cdb->getUrl("TPC_TZERO_OFFSET");
    if (!t0Url.empty())
    {
      CDBTTree tree(t0Url);
      tree.LoadCalibrations();
      result.tzero_ns = tree.GetSingleFloatValue("tpc_tzero");
    }

    return result;
  }

  double TownsendScaledField(const RunConditions &run,
                             const RunConditions &ref)
  {
    if (!run.valid || !ref.valid || run.pressure <= 0.0 || ref.pressure <= 0.0 ||
        run.temperatureK <= 0.0 || ref.temperatureK <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    // Keep E/N fixed: E*T/P = const.
    // E_new = E_ref * (P_new/P_ref) * (T_ref/T_new).
    return kReferenceCMField_Vcm *
           (run.pressure / ref.pressure) *
           (ref.temperatureK / run.temperatureK);
  }

  double T0CorrectedCdbVelocity(const double cdbVelocity_cm_ns,
                                const double cdbT0_ns,
                                const double measuredT0_8bin,
                                const double samplePeriod_ns,
                                const double driftLength_cm)
  {
    if (!std::isfinite(cdbVelocity_cm_ns) || cdbVelocity_cm_ns <= 0.0 ||
        !std::isfinite(cdbT0_ns) || !std::isfinite(measuredT0_8bin) ||
        samplePeriod_ns <= 0.0 || driftLength_cm <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    // User t0 is supplied in units of 8 SAMPA time bins.
    const double measuredT0_ns = measuredT0_8bin * 8.0 * samplePeriod_ns;
    const double deltaT0_ns = measuredT0_ns - cdbT0_ns;

    // Sign convention used here:
    // arrival_time = t0 + drift_time.
    // Therefore a larger t0 means a shorter inferred drift time for the same
    // CM endpoint, and hence a larger corrected drift velocity.
    const double cdbDriftTime_ns = driftLength_cm / cdbVelocity_cm_ns;
    const double correctedDriftTime_ns = cdbDriftTime_ns - deltaT0_ns;

    if (correctedDriftTime_ns <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return driftLength_cm / correctedDriftTime_ns;
  }

  struct SideHistograms
  {
    TH1D *hVzSigned{nullptr};
    TH1D *hVzOut{nullptr};
    TH1D *hVmag{nullptr};
    TH1D *hVr{nullptr};
    TH1D *hVphi{nullptr};
    TProfile *pVzVsR{nullptr};
    TProfile *pVzVsZ{nullptr};
    TProfile *pTrajectoryVzVsR{nullptr};
    TProfile2D *pVzRZ{nullptr};
    TProfile2D *pVzRPhi{nullptr};
  };

  SideHistograms MakeSideHistograms(const int side)
  {
    const char *label = side == 0 ? "south" : "north";
    SideHistograms h;

    h.hVzSigned = new TH1D(Form("h_vz_signed_%s", label),
                           Form("Local v_{z}, %s;v_{z} [cm/#mus];steps", label),
                           500, -10.0, 10.0);
    h.hVzOut = new TH1D(Form("h_vz_out_%s", label),
                        Form("Local outward v_{z}, %s;v_{z,out} [cm/#mus];steps", label),
                        500, 0.0, 10.0);
    h.hVmag = new TH1D(Form("h_vmag_%s", label),
                       Form("Local |v|, %s;|v| [cm/#mus];steps", label),
                       500, 0.0, 10.0);
    h.hVr = new TH1D(Form("h_vr_%s", label),
                     Form("Local v_{r}, %s;v_{r} [cm/#mus];steps", label),
                     500, -2.0, 2.0);
    h.hVphi = new TH1D(Form("h_vphi_%s", label),
                       Form("Local v_{#phi}, %s;v_{#phi} [cm/#mus];steps", label),
                       500, -2.0, 2.0);

    h.pVzVsR = new TProfile(Form("p_vz_vs_r_%s", label),
                            Form("Mean outward v_{z} vs r, %s;r [cm];<v_{z,out}> [cm/#mus]", label),
                            120, 18.0, 82.0);
    h.pVzVsZ = new TProfile(Form("p_vz_vs_z_%s", label),
                            Form("Mean outward v_{z} vs z, %s;z [cm];<v_{z,out}> [cm/#mus]", label),
                            240, -120.0, 120.0);
    h.pTrajectoryVzVsR = new TProfile(Form("p_trajectory_vz_vs_r_%s", label),
                                      Form("CM-to-pad average v_{z} vs pad radius, %s;r_{pad} [cm];<v_{z,out}> [cm/#mus]", label),
                                      120, 18.0, 82.0);

    h.pVzRZ = new TProfile2D(Form("p_vz_rz_%s", label),
                             Form("Mean outward v_{z}, %s;r [cm];z [cm];<v_{z,out}> [cm/#mus]", label),
                             120, 18.0, 82.0,
                             240, -120.0, 120.0);
    h.pVzRPhi = new TProfile2D(Form("p_vz_rphi_%s", label),
                               Form("Mean outward v_{z}, %s;r [cm];#phi [rad];<v_{z,out}> [cm/#mus]", label),
                               120, 18.0, 82.0,
                               144, -TMath::Pi(), TMath::Pi());
    return h;
  }

  bool FillFromPolyline(const TPolyLine3D *poly,
                        const int side,
                        const double step_ns,
                        const double padRadius,
                        SideHistograms &h,
                        double &sumVzOut_cm_us,
                        long long &nSteps)
  {
    if (!poly || poly->GetN() < 2 || !poly->GetP() || step_ns <= 0.0)
    {
      return false;
    }

    const int n = poly->GetN();
    const float *p = poly->GetP();
    const double sideSign = side == 1 ? +1.0 : -1.0;

    double trajectorySumVzOut = 0.0;
    int trajectorySteps = 0;

    for (int i = 0; i < n - 1; ++i)
    {
      const double x0 = p[3 * i + 0];
      const double y0 = p[3 * i + 1];
      const double z0 = p[3 * i + 2];
      const double x1 = p[3 * (i + 1) + 0];
      const double y1 = p[3 * (i + 1) + 1];
      const double z1 = p[3 * (i + 1) + 2];

      // ReverseDrift advances p_{i+1} = p_i - v*dt.
      // Therefore (p_i - p_{i+1})/dt reconstructs the physical electron
      // drift velocity, i.e. CM -> pad plane.
      const double vx = (x0 - x1) / step_ns;  // cm/ns
      const double vy = (y0 - y1) / step_ns;
      const double vz = (z0 - z1) / step_ns;

      const double xm = 0.5 * (x0 + x1);
      const double ym = 0.5 * (y0 + y1);
      const double zm = 0.5 * (z0 + z1);
      const double r = std::hypot(xm, ym);
      const double phi = std::atan2(ym, xm);

      double vr = 0.0;
      double vphi = 0.0;
      if (r > 0.0)
      {
        vr = (vx * xm + vy * ym) / r;
        vphi = (-vx * ym + vy * xm) / r;
      }

      const double vmag = std::sqrt(vx * vx + vy * vy + vz * vz);
      const double vzOut = sideSign * vz;

      // Histograms are in cm/us for readability.
      const double scale = 1000.0;
      h.hVzSigned->Fill(vz * scale);
      h.hVzOut->Fill(vzOut * scale);
      h.hVmag->Fill(vmag * scale);
      h.hVr->Fill(vr * scale);
      h.hVphi->Fill(vphi * scale);
      h.pVzVsR->Fill(r, vzOut * scale);
      h.pVzVsZ->Fill(zm, vzOut * scale);
      h.pVzRZ->Fill(r, zm, vzOut * scale);
      h.pVzRPhi->Fill(r, phi, vzOut * scale);

      trajectorySumVzOut += vzOut;
      ++trajectorySteps;
      sumVzOut_cm_us += vzOut * scale;
      ++nSteps;
    }

    if (trajectorySteps > 0)
    {
      h.pTrajectoryVzVsR->Fill(padRadius,
                               1000.0 * trajectorySumVzOut / trajectorySteps);
    }
    return trajectorySteps > 0;
  }

  void WriteSideHistograms(const SideHistograms &h)
  {
    h.hVzSigned->Write();
    h.hVzOut->Write();
    h.hVmag->Write();
    h.hVr->Write();
    h.hVphi->Write();
    h.pVzVsR->Write();
    h.pVzVsZ->Write();
    h.pTrajectoryVzVsR->Write();
    h.pVzRZ->Write();
    h.pVzRPhi->Write();
  }
}  // namespace TpcDVScan

void Fun4All_TpcDriftVelocityScan(
    int runnumber = 79513,
    const std::string &inputFile =
        "/sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc00_0/run_00079500_00079600/DST_STREAMING_EVENT_ebdc00_0_run3pp_ana532_nocdbtag_v001-00079513-00000.root",
    const double measuredT0_8bin = std::numeric_limits<double>::quiet_NaN(),
    const int nPhi = 24,
    const double reverseDriftStep_ns = 56.8,
    const std::string &outputDir = ".",
    const std::string &globalTag = "newcdbtag")
{
  using namespace TpcDVScan;

  if (inputFile.empty())
  {
    std::cerr << "Need one DST input file from run " << runnumber
              << " containing the RUN geometry nodes (in particular TPCGEOMCONTAINER)."
              << std::endl;
    return;
  }
  if (nPhi < 1 || reverseDriftStep_ns <= 0.0)
  {
    std::cerr << "Invalid nPhi or reverseDriftStep_ns" << std::endl;
    return;
  }

  // -------------------------------------------------------------------------
  // 1. Read run conditions and reconstruction constants directly from CDB.
  // -------------------------------------------------------------------------
  const RunConditions refConditions = ReadAverageRunConditions(kReferenceRun, globalTag);
  const RunConditions runConditions = ReadAverageRunConditions(runnumber, globalTag);
  const CdbRecoConstants cdbReco = ReadCdbRecoConstants(runnumber, globalTag);

  // Restore the requested run after querying the reference run.
  SetCdbRun(runnumber, globalTag);

  double cmField_Vcm = TownsendScaledField(runConditions, refConditions);
  if (!std::isfinite(cmField_Vcm))
  {
    std::cout << "WARNING: could not compute Townsend P/T scaling. Using reference field "
              << kReferenceCMField_Vcm << " V/cm." << std::endl;
    cmField_Vcm = kReferenceCMField_Vcm;
  }

  std::cout << std::setprecision(10)
            << "\n================ run " << runnumber << " ================\n"
            << "reference run = " << kReferenceRun << "\n"
            << "P_ref = " << refConditions.pressure
            << ", T_ref(raw) = " << refConditions.temperatureRaw
            << ", T_ref(K) = " << refConditions.temperatureK << "\n"
            << "P_run = " << runConditions.pressure
            << ", T_run(raw) = " << runConditions.temperatureRaw
            << ", T_run(K) = " << runConditions.temperatureK << "\n"
            << "CM field after Townsend scaling = " << cmField_Vcm << " V/cm\n"
            << "CDB drift velocity = " << cdbReco.driftVelocity_cm_ns
            << " cm/ns = " << 1000.0 * cdbReco.driftVelocity_cm_ns << " cm/us\n"
            << "CDB t0 = " << cdbReco.tzero_ns << " ns\n"
            << "====================================================\n"
            << std::endl;

  // -------------------------------------------------------------------------
  // 2. Minimal Fun4All setup.
  //
  //    - one event DST is used only to drive Fun4All/InitRun
  //    - Tracking_Geometry is loaded separately from CDB into the RUN node,
  //      exactly as in the normal TPC production macro
  //    - PHGarfield is explicitly told NOT to use survey geometry for its
  //      TPC translation/rotation; the loaded TPCGEOMCONTAINER is used only
  //      below to obtain the detector dimensions
  // -------------------------------------------------------------------------
  auto *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // Event input: enough to let Fun4All execute one event and call InitRun.
  auto *in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(inputFile);
  se->registerInputManager(in);

  // Load the standard tracking geometry into the RUN node.
  // This is the same pattern used by the TPC production macro and provides
  // TPCGEOMCONTAINER independently of what is stored in the event DST.
  const std::string geometryUrl =
      CDBInterface::instance()->getUrl("Tracking_Geometry");

  if (geometryUrl.empty())
  {
    std::cerr << "Tracking_Geometry is not available from CDB for run "
              << runnumber << std::endl;
    return;
  }

  std::cout << "Tracking geometry: " << geometryUrl << std::endl;

  auto *geoIn = new Fun4AllRunNodeInputManager("GeoIn");
  geoIn->AddFile(geometryUrl);
  se->registerInputManager(geoIn);

  // Create the standard magnetic-field node without loading the full tracking
  // macro stack.
  const std::string fieldUrl =
      CDBInterface::instance()->getUrl("FIELDMAP_TRACKING");

  if (fieldUrl.empty())
  {
    std::cerr << "FIELDMAP_TRACKING is not available for run "
              << runnumber << std::endl;
    return;
  }

  PHFieldConfigv1 fieldConfig;
  fieldConfig.set_field_config(
      PHFieldConfig::FieldConfigTypes::Field3DCartesian);
  fieldConfig.set_filename(fieldUrl);
  fieldConfig.set_magfield_rescale(1.0);
  PHFieldUtility::GetFieldMapNode(&fieldConfig, se->topNode());

  auto *phg = new PHGarfield("PHGarfield_DVScan");

  // Keep PHGarfield in the nominal/local TPC coordinate system. We still
  // load Tracking_Geometry above so the scan can read the real TPC dimensions,
  // but PHGarfield itself does not apply survey translation/rotation.
  phg->SetUseSurveyGeometry(false);

  phg->SetCMVoltageDefault(cmField_Vcm);

  // This is a run-averaged P/T scan, not an event-BCO scan. Disable the
  // event-current correction explicitly. The base CDB kEff, 3D field maps,
  // frame maps, field-cage setup and gas file are still loaded by PHGarfield.
  phg->SetUseBCOkEffs(false);
  phg->Verbosity(0);
  se->registerSubsystem(phg);

  // One event is sufficient to load the RUN node and execute PHGarfield::InitRun.
  se->run(1);

  // Sanity check that Tracking_Geometry populated the RUN node.
  se->Print("NODETREE");

  // -------------------------------------------------------------------------
  // 3. Get actual TPC z geometry from the run node.
  // -------------------------------------------------------------------------
  auto *tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(se->topNode(), "TPCGEOMCONTAINER");
  if (!tpcGeom)
  {
    std::cerr << "TPCGEOMCONTAINER is missing after loading Tracking_Geometry from CDB." << std::endl;
    se->End();
    delete se;
    return;
  }

  PHG4TpcGeom *layerGeom = tpcGeom->GetLayerCellGeom(20);
  if (!layerGeom)
  {
    std::cerr << "TPC layer-20 geometry is missing." << std::endl;
    se->End();
    delete se;
    return;
  }

  const double maxDriftLength_cm = layerGeom->get_max_driftlength();
  const double cmHalfWidth_cm = layerGeom->get_CM_halfwidth();
  const double zPadNorth_cm = +(maxDriftLength_cm + cmHalfWidth_cm);
  const double zPadSouth_cm = -(maxDriftLength_cm + cmHalfWidth_cm);

  std::cout << std::setprecision(10)
            << "TPC geometry from Tracking_Geometry:"
            << "\n  max drift length = " << maxDriftLength_cm << " cm"
            << "\n  CM half width    = " << cmHalfWidth_cm << " cm"
            << "\n  z pad north      = " << zPadNorth_cm << " cm"
            << "\n  z pad south      = " << zPadSouth_cm << " cm"
            << std::endl;

  const double measuredT0ForCorrection_8bin =
      (std::isfinite(measuredT0_8bin) && measuredT0_8bin > -900.0)
          ? measuredT0_8bin
          : std::numeric_limits<double>::quiet_NaN();

  const double cdbVelocityT0Corrected_cm_ns =
      T0CorrectedCdbVelocity(cdbReco.driftVelocity_cm_ns,
                             cdbReco.tzero_ns,
                             measuredT0ForCorrection_8bin,
                             reverseDriftStep_ns,
                             maxDriftLength_cm);

  if (std::isfinite(measuredT0ForCorrection_8bin))
  {
    const double measuredT0_ns = measuredT0ForCorrection_8bin * 8.0 * reverseDriftStep_ns;
    std::cout << "Measured t0 = " << measuredT0_8bin << " x 8 time bins = "
              << measuredT0_ns << " ns\n"
              << "t0-corrected CDB drift velocity = "
              << cdbVelocityT0Corrected_cm_ns << " cm/ns = "
              << 1000.0 * cdbVelocityT0Corrected_cm_ns << " cm/us"
              << std::endl;
  }

  // -------------------------------------------------------------------------
  // 4. Output histograms and all-step velocity maps.
  // -------------------------------------------------------------------------
  const std::string outputName =
    "/sphenix/tg/tg01/hf/mitrankov/dv/rootfiles/tpc_dv_scan_run" +
    std::to_string(runnumber) + ".root";
  std::unique_ptr<TFile> output(TFile::Open(outputName.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "Could not create " << outputName << std::endl;
    se->End();
    delete se;
    return;
  }

  SideHistograms hSouth = MakeSideHistograms(0);
  SideHistograms hNorth = MakeSideHistograms(1);

  double sumVzOutSouth_cm_us = 0.0;
  double sumVzOutNorth_cm_us = 0.0;
  long long nStepsSouth = 0;
  long long nStepsNorth = 0;
  long long nTrajectoriesSouth = 0;
  long long nTrajectoriesNorth = 0;
  long long nFailedSouth = 0;
  long long nFailedNorth = 0;

  for (int layer = 0; layer < kNLayers; ++layer)
  {
    const double rPad = phg->GetRadius(layer);

    for (int iphi = 0; iphi < nPhi; ++iphi)
    {
      const double phi = 2.0 * TMath::Pi() * (iphi + 0.5) / nPhi - TMath::Pi();
      const double x = rPad * std::cos(phi);
      const double y = rPad * std::sin(phi);

      PHGarfield::ReverseDriftStatus statusNorth;
      std::unique_ptr<TPolyLine3D> north(
          phg->ReverseDrift(x, y, zPadNorth_cm, reverseDriftStep_ns, &statusNorth));
      if (FillFromPolyline(north.get(), 1, reverseDriftStep_ns, rPad,
                           hNorth, sumVzOutNorth_cm_us, nStepsNorth))
      {
        ++nTrajectoriesNorth;
      }
      else
      {
        ++nFailedNorth;
      }

      PHGarfield::ReverseDriftStatus statusSouth;
      std::unique_ptr<TPolyLine3D> south(
          phg->ReverseDrift(x, y, zPadSouth_cm, reverseDriftStep_ns, &statusSouth));
      if (FillFromPolyline(south.get(), 0, reverseDriftStep_ns, rPad,
                           hSouth, sumVzOutSouth_cm_us, nStepsSouth))
      {
        ++nTrajectoriesSouth;
      }
      else
      {
        ++nFailedSouth;
      }
    }
  }

  double garfieldVzSouth_cm_us = nStepsSouth > 0
                                     ? sumVzOutSouth_cm_us / static_cast<double>(nStepsSouth)
                                     : std::numeric_limits<double>::quiet_NaN();

  double garfieldVzNorth_cm_us = nStepsNorth > 0
                                     ? sumVzOutNorth_cm_us / static_cast<double>(nStepsNorth)
                                     : std::numeric_limits<double>::quiet_NaN();

  double garfieldVzCombined_cm_us = (nStepsSouth + nStepsNorth) > 0
                                        ? (sumVzOutSouth_cm_us + sumVzOutNorth_cm_us) /
                                              static_cast<double>(nStepsSouth + nStepsNorth)
                                        : std::numeric_limits<double>::quiet_NaN();

  output->cd();
  WriteSideHistograms(hSouth);
  WriteSideHistograms(hNorth);

  TTree summary("run_summary", "TPC PHGarfield drift-velocity run summary");
  int outRun = runnumber;
  int refRun = kReferenceRun;
  int outNPhi = nPhi;
  int nPTSamples = runConditions.nSamples;
  double pRun = runConditions.pressure;
  double tRunRaw = runConditions.temperatureRaw;
  double tRunK = runConditions.temperatureK;
  double pRef = refConditions.pressure;
  double tRefRaw = refConditions.temperatureRaw;
  double tRefK = refConditions.temperatureK;
  double field = cmField_Vcm;
  double cdbDV = cdbReco.driftVelocity_cm_ns;
  double cdbT0 = cdbReco.tzero_ns;
  double measuredT0 = measuredT0ForCorrection_8bin;
  double correctedDV = cdbVelocityT0Corrected_cm_ns;
  double stepNs = reverseDriftStep_ns;
  double driftLength = maxDriftLength_cm;
  double zPadN = zPadNorth_cm;
  double zPadS = zPadSouth_cm;

  summary.Branch("run", &outRun);
  summary.Branch("reference_run", &refRun);
  summary.Branch("n_phi", &outNPhi);
  summary.Branch("n_pt_samples", &nPTSamples);
  summary.Branch("pressure", &pRun);
  summary.Branch("temperature_raw", &tRunRaw);
  summary.Branch("temperature_K", &tRunK);
  summary.Branch("reference_pressure", &pRef);
  summary.Branch("reference_temperature_raw", &tRefRaw);
  summary.Branch("reference_temperature_K", &tRefK);
  summary.Branch("cm_field_Vcm", &field);
  summary.Branch("cdb_drift_velocity_cm_ns", &cdbDV);
  summary.Branch("cdb_tzero_ns", &cdbT0);
  summary.Branch("measured_t0_8bin", &measuredT0);
  summary.Branch("cdb_drift_velocity_t0corr_cm_ns", &correctedDV);
  summary.Branch("reverse_drift_step_ns", &stepNs);
  summary.Branch("max_drift_length_cm", &driftLength);
  summary.Branch("z_pad_north_cm", &zPadN);
  summary.Branch("z_pad_south_cm", &zPadS);
  summary.Branch("garfield_mean_vz_south_cm_us", &garfieldVzSouth_cm_us);
  summary.Branch("garfield_mean_vz_north_cm_us", &garfieldVzNorth_cm_us);
  summary.Branch("garfield_mean_vz_combined_cm_us", &garfieldVzCombined_cm_us);
  summary.Branch("n_steps_south", &nStepsSouth);
  summary.Branch("n_steps_north", &nStepsNorth);
  summary.Branch("n_trajectories_south", &nTrajectoriesSouth);
  summary.Branch("n_trajectories_north", &nTrajectoriesNorth);
  summary.Branch("n_failed_south", &nFailedSouth);
  summary.Branch("n_failed_north", &nFailedNorth);
  summary.Fill();
  summary.Write();

  TParameter<double>("reference_cm_field_Vcm", kReferenceCMField_Vcm).Write();
  TParameter<int>("reference_run", kReferenceRun).Write();

  output->Write();
  output->Close();

  std::cout << std::setprecision(8)
            << "\nPHGarfield all-step mean outward vz:\n"
            << "  South   = " << garfieldVzSouth_cm_us << " cm/us\n"
            << "  North   = " << garfieldVzNorth_cm_us << " cm/us\n"
            << "  Combined= " << garfieldVzCombined_cm_us << " cm/us\n"
            << "CDB       = " << 1000.0 * cdbReco.driftVelocity_cm_ns << " cm/us\n";
  if (std::isfinite(cdbVelocityT0Corrected_cm_ns))
  {
    std::cout << "CDB(t0 corrected) = "
              << 1000.0 * cdbVelocityT0Corrected_cm_ns << " cm/us\n";
  }
  std::cout << "Saved: " << outputName << std::endl;

  se->End();
  delete se;
}
