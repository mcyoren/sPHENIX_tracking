// TPC-only raw-hit unpacking and PHGarfield drift-velocity QA.
// Keep GlobalVariables.C first, with a blank line after it.
#include <GlobalVariables.C>

#include <G4_TrkrVariables.C>
#include <Trkr_TpcReadoutInit.C>

#include <phgarfielddriftvelocityqa/PHGarfieldDriftVelocityQA.h>
#include <qautils/QAHistManagerDef.h>
#include <tpc/TpcCombinedRawDataUnpacker.h>
#include <tpcconditions/TpcConditionsReco.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/recoConsts.h>

#include <TSystem.h>

#include <format>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libTpcConditions.so)
R__LOAD_LIBRARY(libPHGarfieldDriftVelocityQA.so)

namespace DriftVelocityQA
{
  class SkipFirstN : public SubsysReco
  {
   public:
    explicit SkipFirstN(const int n)
      : SubsysReco("DriftVelocityQA_SkipFirstN")
      , m_target(n)
    {
    }

    int process_event(PHCompositeNode* /*topNode*/) override
    {
      if (m_count < m_target)
      {
        ++m_count;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }

   private:
    int m_target = 0;
    int m_count = 0;
  };
}  // namespace DriftVelocityQA

void Fun4All_raw_hit_TPC_DriftVelocityQA(
    const int nEvents = 2,
    const int runnumber = 79513,
    const int segment = 0,
    const std::string& outdir = ".",
    const int nSkip = 0,
    const std::string& collision = "run3pp",
    const std::string& production = "ana532_nocdbtag_v001",
    const std::string& outfilename = "HITS_ppFieldOn",
    const std::string& datatype = "physics")
{
  if (nEvents < 0 || nSkip < 0 || runnumber <= 0 || segment < 0)
  {
    std::cerr << "Invalid event count, skip count, run number, or segment" << std::endl;
    gSystem->Exit(1);
    return;
  }

  auto* se = Fun4AllServer::instance();
  se->Verbosity(1);

  auto* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  Enable::CDB = true;

  // Direct histogram registration does not need Enable::QA or QA.C.
  QAHistManagerDef::getHistoManager();
  se->registerSubsystem(new DriftVelocityQA::SkipFirstN(nSkip));
  se->registerSubsystem(new FlagHandler());

  const int runbase = (runnumber / 100) * 100;
  const int runnext = runbase + 100;

  // Read the 24 EBDCs x 2 endpoints, without other detector streams.
  for (int ebdc = 0; ebdc < 24; ++ebdc)
  {
    for (int endpoint = 0; endpoint < 2; ++endpoint)
    {
      const std::string stream = std::format("ebdc{:02d}_{}", ebdc, endpoint);
      const std::string filename = std::format(
          "DST_STREAMING_EVENT_{}_{}_{}-{:08d}-{:05d}.root",
          stream, collision, production, runnumber, segment);
      const std::string filepath = std::format(
          "/sphenix/lustre01/sphnxpro/production/{}/{}/{}/"
          "DST_STREAMING_EVENT_{}/run_{:08d}_{:08d}/{}",
          collision, datatype, production, stream, runbase, runnext, filename);
      const std::string input = gSystem->AccessPathName(filename.c_str()) ? filepath : filename;

      std::cout << "Adding DST: " << input << std::endl;
      auto* hitsin = new Fun4AllDstInputManager("InputManager_" + stream);
      if (hitsin->fileopen(input) != 0)
      {
        std::cerr << "Failed to open " << input << std::endl;
        delete hitsin;
        delete se;
        gSystem->Exit(1);
        return;
      }
      se->registerInputManager(hitsin);
    }
  }

  // Geometry is supplied by CDB; no ACTS geometry or tracking initialization.
  const std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  if (geofile.empty())
  {
    std::cerr << "Missing Tracking_Geometry CDB payload" << std::endl;
    delete se;
    gSystem->Exit(1);
    return;
  }
  auto* ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TRACKING::streaming_mode = true;
  TRACKING::tpc_zero_supp = true;
  TpcReadoutInit(runnumber);

  std::cout << "Run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " presample shift: " << TRACKING::reco_tpc_time_presample
            << " unpacker t0: " << TRACKING::reco_t0 << std::endl;

  for (int ebdc = 0; ebdc < 24; ++ebdc)
  {
    for (int endpoint = 0; endpoint < 2; ++endpoint)
    {
      const std::string suffix = std::format("{:02d}_{}", ebdc, endpoint);
      auto* unpacker = new TpcCombinedRawDataUnpacker("TpcCombinedRawDataUnpacker" + suffix);
      unpacker->set_presampleShift(TRACKING::reco_tpc_time_presample);
      unpacker->set_t0(TRACKING::reco_t0);
      unpacker->useRawHitNodeName("TPCRAWHIT_" + suffix);
      unpacker->ReadZeroSuppressedData();
      unpacker->Verbosity(0);
      se->registerSubsystem(unpacker);
    }
  }

  se->registerSubsystem(new TpcConditionsReco());

  auto* qa = new PHGarfieldDriftVelocityQA();
  qa->setUsePTCorrectedPHGarfield(true);
  qa->setReferenceCMFieldVcm(377.0);
  qa->setReferencePressure(30.0968906);
  qa->setReferenceTemperatureK(299.98646566);
  qa->setTpcTimeBinNs(56.8);
  qa->setPHGarfieldT0Bins(8.0);
  qa->setPHGarfieldStepNs(5.0);
  se->registerSubsystem(qa);

  const std::string qadir = outdir + "/QA";
  if (gSystem->AccessPathName(qadir.c_str()) && gSystem->mkdir(qadir.c_str(), true) != 0)
  {
    std::cerr << "Cannot create QA directory: " << qadir << std::endl;
    delete se;
    gSystem->Exit(1);
    return;
  }

  // Fun4All run(0) processes all available events, including after nSkip.
  const int status = se->run(nEvents == 0 ? 0 : nEvents + nSkip);
  const int endStatus = se->End();

  const std::string output = std::format(
      "{}/tpc_drift_velocity_{}_{}evt_{}skip_{}_{}_qa.root",
      qadir, outfilename, nEvents, nSkip, runnumber, segment);
  QAHistManagerDef::saveQARootFile(output);
  std::cout << "QA output: " << output << std::endl;

  se->PrintTimer();
  CDBInterface::instance()->Print();
  delete se;
  gSystem->Exit(status == 0 && endStatus == 0 ? 0 : 1);
}
