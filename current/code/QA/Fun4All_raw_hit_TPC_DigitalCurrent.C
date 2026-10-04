#include <fun4allraw/TpcTimeFrameBuilderRun3.h>
#include <qautils/QAHistManagerDef.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <phool/recoConsts.h>

#include <Rtypes.h>
#include <TSystem.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libEvent.so)

namespace
{
  std::string raw_file_name(
      const std::string& rawdir,
      const int ebdc,
      const int endpoint,
      const int runnumber,
      const int segment)
  {
    std::ostringstream s;
    s << rawdir;
    if (!rawdir.empty() && rawdir.back() != '/') s << "/";
    s << "TPC_ebdc"
      << std::setw(2) << std::setfill('0') << ebdc
      << "_" << endpoint
      << "_physics-"
      << std::setw(8) << std::setfill('0') << runnumber
      << "-"
      << std::setw(4) << std::setfill('0') << segment
      << ".evt";
    return s.str();
  }

  bool is_tpc_packet(const int packetid)
  {
    if (packetid < 4000 || packetid > 4231) return false;
    const int endpoint = packetid % 10;
    const int ebdc = (packetid - 4000) / 10;
    return ebdc >= 0 && ebdc < 24 && (endpoint == 0 || endpoint == 1);
  }

  std::string output_name(
      const std::string& outdir,
      const std::string& outfilename,
      const int runnumber,
      const int segment,
      const int packetid)
  {
    std::ostringstream s;
    s << outdir << "/" << outfilename
      << "_DigitalCurrent_"
      << std::setw(8) << std::setfill('0') << runnumber << "_"
      << std::setw(4) << std::setfill('0') << segment
      << "_packet" << packetid << ".root";
    return s.str();
  }
}

/*
 * Extract firmware TPC Digital Current directly from modern raw endpoint files.
 *
 * First 8 arguments match Fun4All_raw_hit_TPC_reco_Work.C.
 *
 * Raw files are found automatically as
 *
 * /sphenix/lustre01/sphnxpro/fromhpss/physics_2025/tpc/physics/
 * TPC_ebdcXX_Y_physics-000RUNNN-SSSS.evt
 *
 * There are at most 50 DATAEVENTs per endpoint raw file in this workflow,
 * so the macro processes min(nEvents, 50-nSkip) DATAEVENTs per endpoint.
 *
 * Example:
 *
 * .x Fun4All_raw_hit_TPC_DigitalCurrent.C(
 *      10,81558,0,
 *      "/sphenix/tg/tg01/hf/mitrankov/ibf/output/",
 *      0,"run3pp","ana532_nocdbtag_v001","RawHitQA")
 */
int Fun4All_raw_hit_TPC_DigitalCurrent(
    const int nEvents = 10,
    const int runnumber = 81558,
    const int segment = 0,
    const std::string& outdir = ".",
    const int nSkip = 0,
    const std::string& collision = "run3pp",
    const std::string& production = "ana532_nocdbtag_v001",
    const std::string& outfilename = "RawHitQA",
    const std::string& rawdir =
        "/sphenix/lustre01/sphnxpro/fromhpss/physics_2025/tpc/physics",
    const int verbosity = 0)
{
  (void) collision;
  (void) production;

  constexpr int MAX_EVENTS_PER_RAW_ENDPOINT = 50;
  const int available_after_skip = std::max(0, MAX_EVENTS_PER_RAW_ENDPOINT - nSkip);
  const int events_per_endpoint =
      (nEvents < 0) ? available_after_skip : std::min(nEvents, available_after_skip);

  if (events_per_endpoint <= 0)
  {
    std::cout << "Nothing to process: nSkip=" << nSkip
              << ", maximum DATAEVENTs per endpoint="
              << MAX_EVENTS_PER_RAW_ENDPOINT << std::endl;
    return 0;
  }

  std::cout << "Requested nEvents = " << nEvents
            << ", nSkip = " << nSkip
            << " -> processing up to " << events_per_endpoint
            << " DATAEVENTs per endpoint" << std::endl;

  auto* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  QAHistManagerDef::getHistoManager();

  std::string dcdir = outdir;
  if (!dcdir.empty() && dcdir.back() == '/') dcdir.pop_back();
  dcdir += "/DC";
  gSystem->mkdir(dcdir.c_str(), true);

  std::map<int, std::unique_ptr<TpcTimeFrameBuilderRun3>> builders;
  std::map<int, long long> packet_calls;

  int files_found = 0;
  int files_missing = 0;
  long long data_events_seen = 0;
  long long data_events_used = 0;

  // 24 EBDCs x 2 endpoints = 48 raw files for one raw segment.
  for (int ebdc = 0; ebdc < 24; ++ebdc)
  {
    for (int endpoint = 0; endpoint < 2; ++endpoint)
    {
      const std::string rawfile =
          raw_file_name(rawdir, ebdc, endpoint, runnumber, segment);

      if (gSystem->AccessPathName(rawfile.c_str()))
      {
        ++files_missing;
        std::cout << "MISSING: " << rawfile << std::endl;
        continue;
      }

      ++files_found;
      std::cout << "\nEBDC "
                << std::setw(2) << std::setfill('0') << ebdc
                << " endpoint " << endpoint
                << "\n  " << rawfile << std::endl;

      int status = 0;
      std::unique_ptr<Eventiterator> it(
          new fileEventiterator(rawfile.c_str(), status));

      if (status != 0 || !it)
      {
        std::cout << "  ERROR: fileEventiterator status=" << status << std::endl;
        continue;
      }

      int skipped = 0;
      int used = 0;

      while (used < events_per_endpoint)
      {
        Event* evt = it->getNextEvent();
        if (!evt) break;
        std::unique_ptr<Event> event(evt);

        if (event->getEvtType() != DATAEVENT) continue;
        ++data_events_seen;

        if (skipped < nSkip)
        {
          ++skipped;
          continue;
        }

        constexpr int MAX_PACKETS = 128;
        Packet* packets[MAX_PACKETS] = {nullptr};
        const int npackets = event->getPacketList(packets, MAX_PACKETS);

        for (int ip = 0; ip < npackets; ++ip)
        {
          Packet* packet = packets[ip];
          if (!packet) continue;

          const int packetid = packet->getIdentifier();

          if (is_tpc_packet(packetid))
          {
            auto found = builders.find(packetid);

            if (found == builders.end())
            {
              auto builder =
                  std::make_unique<TpcTimeFrameBuilderRun3>(packetid);
              builder->setVerbosity(verbosity);

              const std::string outfile =
                  output_name(dcdir, outfilename, runnumber, segment, packetid);

              builder->SaveDigitalCurrentDebugTTree(outfile);

              std::cout << "  discovered packet " << packetid
                        << "\n    -> " << outfile << std::endl;

              found = builders.emplace(packetid, std::move(builder)).first;
            }

            found->second->ProcessPacket(packet);
            ++packet_calls[packetid];
          }
          else if (verbosity > 1)
          {
            std::cout << "  ignoring packet " << packetid << std::endl;
          }

          delete packet;
          packets[ip] = nullptr;
        }

        ++used;
        ++data_events_used;
      }

      std::cout << "  processed " << used
                << " DATAEVENTs after skipping " << skipped << std::endl;
    }
  }

  std::cout << "\n============================================================\n";
  std::cout << "Digital Current extraction summary\n";
  std::cout << "  raw endpoint files found : " << files_found << "/48\n";
  std::cout << "  raw endpoint files missing: " << files_missing << "/48\n";
  std::cout << "  DATAEVENTs seen           : " << data_events_seen << "\n";
  std::cout << "  DATAEVENTs processed      : " << data_events_used << "\n";
  std::cout << "  TPC packet IDs seen       : " << builders.size() << "\n";

  for (const auto& [packetid, count] : packet_calls)
  {
    std::cout << "    packet " << packetid
              << ": ProcessPacket calls = " << count << "\n";
  }
  std::cout << "============================================================\n";

  if (files_found == 0)
  {
    std::cout
        << "\nERROR: no raw files found for run " << runnumber
        << ", raw segment " << segment << ".\n"
        << "Expected names like:\n  "
        << raw_file_name(rawdir, 0, 0, runnumber, segment) << "\n";
    return 2;
  }

  if (builders.empty())
  {
    std::cout
        << "\nERROR: raw files were opened but no TPC packet IDs were decoded.\n"
        << "Run once with the optional final argument verbosity=2 if needed.\n";
    return 3;
  }

  // T_DigitalCurrent is written/flushed by the builder destructor.
  builders.clear();

  std::cout << "\nFinished. Output ROOT files:\n  "
            << dcdir << "\n"
            << "Each packet file should contain T_DigitalCurrent.\n";

  return 0;
}
