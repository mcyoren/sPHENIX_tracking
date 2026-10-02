#include "T_DigitalCurrent.C"
#include "InfoTPC.h"
#include "DCGemCurrent.h"

R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dcmaps
{
  constexpr ULong64_t FRAME_BCO = 9720ULL;
  constexpr int N_SIDES = 2;
  constexpr int N_SECTORS = 12;
  constexpr int N_MODULES = 3;
  constexpr int N_VARIANTS = 4;
  constexpr int N_FRAME_BINS = 80;
  constexpr double Z_MAX_MM = 1075.0; // follows the existing T_DigitalCurrent.C convention
  constexpr int N_PADS[N_MODULES] = {94, 128, 192}; // matches the analysis notebook
  constexpr double FIXED_PEDESTAL = 60.0;
  constexpr double STRICT_MAX_PEDSUB = 30000.0;

  enum Variant
  {
    RAW = 0,
    PED60_SIGNED = 1,
    PED60_POSITIVE = 2,
    PED60_STRICT = 3
  };

  const char* VARIANT_NAME[N_VARIANTS] = {
      "raw", "ped60_signed", "ped60_positive", "ped60_strict"};

  struct Interval
  {
    ULong64_t begin = 0; // inclusive
    ULong64_t end = 0;   // exclusive
  };

  struct FileRecord
  {
    int run = 0;
    int ebdc = -1;
    int stream = -1;
    int segment = -1;
    ULong64_t first_bco = 0;
    ULong64_t last_bco = 0;
    Long64_t entries = 0;
    std::string path;
  };

  struct PadInfo
  {
    bool valid = false;
    double r_mm = 0.0;
    double phi_base = 0.0;
    int module = -1;
    int local_layer = -999; // negative values are antenna-pad rings in R1
    int local_pad = -1;
  };

  using PadLookup = std::vector<std::vector<std::optional<PadInfo>>>;

  std::vector<std::string> split_csv(const std::string& line)
  {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i)
    {
      const char c = line[i];
      if (c == '"')
      {
        if (in_quotes && i + 1 < line.size() && line[i + 1] == '"')
        {
          cur.push_back('"');
          ++i;
        }
        else
        {
          in_quotes = !in_quotes;
        }
      }
      else if (c == ',' && !in_quotes)
      {
        out.push_back(cur);
        cur.clear();
      }
      else
      {
        cur.push_back(c);
      }
    }
    out.push_back(cur);
    return out;
  }

  int csv_to_input_fee(int module, int csv_idx)
  {
    static const std::vector<int> R1 = {2,4,3,13,17,16};
    static const std::vector<int> R2 = {11,12,19,18,0,1,15,14};
    static const std::vector<int> R3 = {20,22,21,23,25,24,10,9,8,6,7,5};
    const std::vector<int>* v = nullptr;
    if (module == 0) v = &R1;
    if (module == 1) v = &R2;
    if (module == 2) v = &R3;
    if (!v || csv_idx < 0 || csv_idx >= static_cast<int>(v->size())) return -1;
    return (*v)[csv_idx];
  }

  PadLookup build_detailed_pad_lookup(const std::string& mapping_dir)
  {
    PadLookup lookup(26);
    const std::array<std::string,3> names = {
        "R1_ChannelMapping.csv", "R2_ChannelMapping.csv", "R3_ChannelMapping.csv"};

    for (int module = 0; module < N_MODULES; ++module)
    {
      const std::string path = mapping_dir + "/" + names[module];
      std::ifstream in(path);
      if (!in) throw std::runtime_error("Cannot open mapping CSV: " + path);

      std::string header;
      std::getline(in, header);
      const auto h = split_csv(header);

      auto col = [&](const std::string& name) -> int
      {
        for (int i = 0; i < static_cast<int>(h.size()); ++i)
          if (h[i] == name) return i;
        return -1;
      };

      const int cRadius = col("Radius");
      const int cPad = col("Pad");
      const int cFee = col("FEE");
      const int cChan = col("FEE_Chan");
      const int cPadR = col("PadR");
      const int cPadPhi = col("PadPhi");
      if (cRadius < 0 || cPad < 0 || cFee < 0 || cChan < 0 || cPadR < 0 || cPadPhi < 0)
        throw std::runtime_error("Required columns missing from " + path);

      std::string line;
      while (std::getline(in, line))
      {
        if (line.empty()) continue;
        const auto f = split_csv(line);
        const int max_col = std::max({cRadius,cPad,cFee,cChan,cPadR,cPadPhi});
        if (static_cast<int>(f.size()) <= max_col) continue;

        const int csv_fee = static_cast<int>(std::llround(std::stod(f[cFee])));
        const int input_fee = csv_to_input_fee(module, csv_fee);
        if (input_fee < 0 || input_fee >= 26) continue;

        const int channel = std::stoi(f[cChan]);
        if (channel < 0) continue;
        if (static_cast<int>(lookup[input_fee].size()) <= channel)
          lookup[input_fee].resize(channel + 1);

        PadInfo p;
        p.valid = true;
        p.r_mm = std::stod(f[cPadR]);
        p.phi_base = std::stod(f[cPadPhi]);
        p.module = module;
        p.local_layer = static_cast<int>(std::llround(std::stod(f[cRadius])));
        p.local_pad = static_cast<int>(std::llround(std::stod(f[cPad])));
        lookup[input_fee][channel] = p;
      }
    }
    return lookup;
  }

  std::optional<PadInfo> get_pad(const PadLookup& lookup, int fee, int channel)
  {
    if (fee < 0 || fee >= static_cast<int>(lookup.size())) return std::nullopt;
    if (channel < 0 || channel >= static_cast<int>(lookup[fee].size())) return std::nullopt;
    return lookup[fee][channel];
  }

  std::vector<FileRecord> read_index(const std::string& filename, int run)
  {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open index: " + filename);
    std::vector<FileRecord> out;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line))
    {
      if (line.empty()) continue;
      const auto f = split_csv(line);
      if (f.size() < 8) continue;
      FileRecord r;
      r.run = std::stoi(f[0]);
      if (r.run != run) continue;
      r.ebdc = std::stoi(f[1]);
      r.stream = std::stoi(f[2]);
      r.segment = std::stoi(f[3]);
      r.first_bco = std::stoull(f[4]);
      r.last_bco = std::stoull(f[5]);
      r.entries = std::stoll(f[6]);
      r.path = f[7];
      out.push_back(r);
    }
    return out;
  }

  std::vector<Interval> read_window(const std::string& filename, int run, int window)
  {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open windows file: " + filename);
    std::vector<Interval> out;
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line))
    {
      if (line.empty()) continue;
      const auto f = split_csv(line);
      if (f.size() < 6) continue;
      if (std::stoi(f[0]) != run || std::stoi(f[1]) != window) continue;
      out.push_back({std::stoull(f[3]), std::stoull(f[4])});
    }
    std::sort(out.begin(), out.end(), [](const Interval& a, const Interval& b){ return a.begin < b.begin; });
    return out;
  }

  bool overlap(const Interval& a, const Interval& b)
  {
    return std::min(a.end, b.end) > std::max(a.begin, b.begin);
  }

  Interval intersection(const Interval& a, const Interval& b)
  {
    return {std::max(a.begin,b.begin), std::min(a.end,b.end)};
  }

  double wrap_phi(double phi)
  {
    const double twopi = 2.0 * M_PI;
    phi = std::fmod(phi, twopi);
    if (phi < 0) phi += twopi;
    return phi;
  }

  ULong64_t live_midpoint_bco(const std::vector<Interval>& intervals)
  {
    ULong64_t total = 0;
    for (const auto& x : intervals) total += x.end - x.begin;
    if (total == 0 || intervals.empty()) return 0;

    ULong64_t offset = total / 2;
    for (const auto& x : intervals)
    {
      const ULong64_t length = x.end - x.begin;
      if (offset < length) return x.begin + offset;
      offset -= length;
    }
    return intervals.back().end - 1;
  }

  Long64_t lower_bound_bco(T_DigitalCurrent& t, ULong64_t target)
  {
    Long64_t lo = 0;
    Long64_t hi = t.GetEntries();
    while (lo < hi)
    {
      const Long64_t mid = lo + (hi - lo) / 2;
      if (t.LoadTree(mid) < 0) { hi = mid; continue; }
      t.GetEntry(mid);
      const ULong64_t bco = t.dc_gtm_bco;
      if (bco == 0 || bco < target) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  struct Histograms
  {
    std::array<std::array<TH3D*,N_SIDES>,N_VARIANTS> frame{};
    std::array<std::array<TH3D*,N_SIDES>,N_VARIANTS> z{};
    std::array<TH3D*,N_SIDES> exposure_frame{};
    std::array<TH3D*,N_SIDES> count_frame{};

    TH2D* native[N_VARIANTS][N_SIDES][N_SECTORS][N_MODULES]{};
    TH2D* native_exposure[N_SIDES][N_SECTORS][N_MODULES]{};
    TH2D* native_count[N_SIDES][N_SECTORS][N_MODULES]{};
  };

  Histograms make_histograms(int run, int window)
  {
    Histograms H;

    for (int v = 0; v < N_VARIANTS; ++v)
    {
      for (int side = 0; side < N_SIDES; ++side)
      {
        H.frame[v][side] = new TH3D(
            Form("h_dc_phi_r_frame_%s_side%d", VARIANT_NAME[v], side),
            Form("DC %s run %d window %d side %d;#phi;r [mm];DC frame mod 80",
                 VARIANT_NAME[v], run, window, side),
            TPC::nphi, TPC::phi_bins,
            TPC::nr, TPC::r_bins,
            N_FRAME_BINS, -0.5, 79.5);
        H.frame[v][side]->SetDirectory(nullptr);

        const double zlo = side == 0 ? -Z_MAX_MM : 0.0;
        const double zhi = side == 0 ? 0.0 : Z_MAX_MM;
        H.z[v][side] = new TH3D(
            Form("h_dc_phi_r_z_%s_side%d", VARIANT_NAME[v], side),
            Form("DC %s run %d window %d side %d;#phi;r [mm];z [mm]",
                 VARIANT_NAME[v], run, window, side),
            TPC::nphi, TPC::phi_bins,
            TPC::nr, TPC::r_bins,
            N_FRAME_BINS, zlo, zhi);
        H.z[v][side]->SetDirectory(nullptr);
      }
    }

    for (int side = 0; side < N_SIDES; ++side)
    {
      H.exposure_frame[side] = new TH3D(
          Form("h_dc_nsamples_phi_r_frame_side%d", side),
          Form("DC nsamples exposure run %d window %d side %d;#phi;r [mm];DC frame mod 80",
               run, window, side),
          TPC::nphi, TPC::phi_bins,
          TPC::nr, TPC::r_bins,
          N_FRAME_BINS, -0.5, 79.5);
      H.exposure_frame[side]->SetDirectory(nullptr);

      H.count_frame[side] = new TH3D(
          Form("h_dc_sampleblocks_phi_r_frame_side%d", side),
          Form("DC sample-block count run %d window %d side %d;#phi;r [mm];DC frame mod 80",
               run, window, side),
          TPC::nphi, TPC::phi_bins,
          TPC::nr, TPC::r_bins,
          N_FRAME_BINS, -0.5, 79.5);
      H.count_frame[side]->SetDirectory(nullptr);

      for (int sec = 0; sec < N_SECTORS; ++sec)
      {
        for (int mod = 0; mod < N_MODULES; ++mod)
        {
          for (int v = 0; v < N_VARIANTS; ++v)
          {
            H.native[v][side][sec][mod] = new TH2D(
                Form("h_dc_pad_vs_layer_%s_side%d_sec%02d_R%d",
                     VARIANT_NAME[v], side, sec, mod+1),
                Form("DC %s side %d sector %02d R%d;local pad;local layer",
                     VARIANT_NAME[v], side, sec, mod+1),
                N_PADS[mod], -0.5, N_PADS[mod]-0.5,
                16, -0.5, 15.5);
            H.native[v][side][sec][mod]->SetDirectory(nullptr);
          }

          H.native_exposure[side][sec][mod] = new TH2D(
              Form("h_dc_nsamples_pad_vs_layer_side%d_sec%02d_R%d", side, sec, mod+1),
              Form("DC nsamples exposure side %d sector %02d R%d;local pad;local layer",
                   side, sec, mod+1),
              N_PADS[mod], -0.5, N_PADS[mod]-0.5,
              16, -0.5, 15.5);
          H.native_exposure[side][sec][mod]->SetDirectory(nullptr);

          H.native_count[side][sec][mod] = new TH2D(
              Form("h_dc_sampleblocks_pad_vs_layer_side%d_sec%02d_R%d", side, sec, mod+1),
              Form("DC sample-block count side %d sector %02d R%d;local pad;local layer",
                   side, sec, mod+1),
              N_PADS[mod], -0.5, N_PADS[mod]-0.5,
              16, -0.5, 15.5);
          H.native_count[side][sec][mod]->SetDirectory(nullptr);
        }
      }
    }

    return H;
  }

  void fill_variant(Histograms& H, int variant, int side, int sector,
                    const PadInfo& p, int frame, double z_mm, double weight)
  {
    const double phi = wrap_phi(p.phi_base + (side == 0 ? sector - 12 : sector) * M_PI / 6.0);
    H.frame[variant][side]->Fill(phi, p.r_mm, frame, weight);
    H.z[variant][side]->Fill(phi, p.r_mm, z_mm, weight);

    if (p.local_layer >= 0 && p.local_layer < 16 &&
        p.module >= 0 && p.module < N_MODULES &&
        p.local_pad >= 0 && p.local_pad < N_PADS[p.module])
    {
      H.native[variant][side][sector][p.module]->Fill(p.local_pad, p.local_layer, weight);
    }
  }
}

void make_dc_maps(
    int run = 81558,
    int window = 0,
    const char* indexCsv = "dc_index_run81558.csv",
    const char* windowsCsv = "dc_windows_run81558.csv",
    const char* outputRoot = "dc_run81558_win00.root",
    const char* mappingDir = ".",
    const char* cdbGlobalTag = "newcdbtag")
{
  using namespace dcmaps;

  std::vector<FileRecord> files;
  std::vector<Interval> win;
  PadLookup pad_lookup;
  try
  {
    files = read_index(indexCsv, run);
    win = read_window(windowsCsv, run, window);
    pad_lookup = build_detailed_pad_lookup(mappingDir);
  }
  catch (const std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return;
  }

  if (win.empty())
  {
    std::cerr << "ERROR: no intervals for run " << run << " window " << window << std::endl;
    return;
  }

  ULong64_t live_bco_ticks = 0;
  ULong64_t live_frame_slots = 0;
  ULong64_t bco_start = win.front().begin;
  ULong64_t bco_end = win.back().end;
  for (const auto& x : win)
  {
    live_bco_ticks += x.end - x.begin;
    live_frame_slots += (x.end - x.begin) / FRAME_BCO;
  }

  std::cout << "Run " << run << " window " << window
            << ": " << win.size() << " synchronized subinterval(s), "
            << live_frame_slots << " live frame slots" << std::endl;

  // One representative GEM-current value per synchronized output file.
  // Use the midpoint in accumulated live BCO exposure, not simply
  // (bco_start+bco_end)/2, so windows containing gaps are handled correctly.
  const ULong64_t conditions_bco = live_midpoint_bco(win);
  dcgem::Result gem;
  try
  {
    gem = dcgem::load_r1_currents(run, conditions_bco, cdbGlobalTag);
  }
  catch (const std::exception& e)
  {
    std::cerr << "WARNING: GEM current lookup failed: " << e.what() << std::endl;
  }

  Histograms H = make_histograms(run, window);

  // side x metric
  // Last four bins carry hadd-safe GEM-current accumulators.
  std::array<std::array<long double,16>,N_SIDES> norm{};
  std::array<std::array<long double,6>,N_SIDES> rejection{};
  std::array<std::set<ULong64_t>,N_SIDES> observed_slots;

  // Additive quantities that remain valid after hadd.
  for (int side = 0; side < N_SIDES; ++side)
  {
    norm[side][0] = 1.0; // N_WINDOWS
    norm[side][1] = static_cast<long double>(live_bco_ticks);
    norm[side][2] = static_cast<long double>(live_frame_slots);
    norm[side][3] = static_cast<long double>(live_frame_slots) / 80.0L;

    if (gem.available)
    {
      const long double current = side == 0
          ? static_cast<long double>(gem.load_sr1)
          : static_cast<long double>(gem.load_nr1);
      norm[side][12] = current;
      norm[side][13] = current * current;
      norm[side][14] = current * static_cast<long double>(live_frame_slots);
      norm[side][15] = 1.0L;
    }
  }

  // Group index records by EBDC/stream and sort by BCO so overlapping files are clipped.
  std::array<std::vector<FileRecord>,48> by_stream;
  for (const auto& f : files)
  {
    if (f.ebdc < 0 || f.ebdc >= 24 || f.stream < 0 || f.stream >= 2) continue;
    const Interval fi{f.first_bco, f.last_bco + FRAME_BCO};
    bool useful = false;
    for (const auto& w : win) if (overlap(fi,w)) { useful = true; break; }
    if (useful) by_stream[2*f.ebdc + f.stream].push_back(f);
  }

  for (auto& v : by_stream)
  {
    std::sort(v.begin(), v.end(), [](const FileRecord& a, const FileRecord& b)
    {
      if (a.first_bco != b.first_bco) return a.first_bco < b.first_bco;
      return a.segment < b.segment;
    });
  }

  for (int key = 0; key < 48; ++key)
  {
    const int ebdc = key / 2;
    const int stream = key % 2;
    const int side = ebdc >= 12 ? 1 : 0;
    const int sector = ebdc % 12;

    if (by_stream[key].empty())
    {
      std::cerr << "ERROR: no files overlap window for ebdc=" << ebdc
                << " stream=" << stream << std::endl;
      continue;
    }

    ULong64_t processed_until = 0;
    ULong64_t covered_frames = 0;

    for (const auto& f : by_stream[key])
    {
      const Interval file_interval{f.first_bco, f.last_bco + FRAME_BCO};
      std::unique_ptr<T_DigitalCurrent> t(new T_DigitalCurrent(f.path.c_str(), ebdc));

      for (const auto& w : win)
      {
        if (!overlap(file_interval, w)) continue;
        Interval use = intersection(file_interval, w);
        if (processed_until > use.begin) use.begin = processed_until;
        if (use.end <= use.begin) continue;

        const ULong64_t nframes_piece = (use.end - use.begin) / FRAME_BCO;
        covered_frames += nframes_piece;
        processed_until = std::max(processed_until, use.end);

        std::cout << "Processing ebdc=" << ebdc
                  << " stream=" << stream
                  << " seg=" << f.segment
                  << " BCO [" << use.begin << "," << use.end << ")"
                  << " file=" << f.path << std::endl;

        Long64_t first_entry = lower_bound_bco(*t, use.begin);
        const Long64_t nentries = t->GetEntries();

        for (Long64_t j = first_entry; j < nentries; ++j)
        {
          if (t->LoadTree(j) < 0) break;
          t->GetEntry(j);
          const ULong64_t bco = t->dc_gtm_bco;

          if (bco == 0 || bco < use.begin) continue;
          if (bco >= use.end) break;

          if (t->dc_data_crc != t->dc_calc_crc)
          {
            rejection[side][0] += 1.0L; // BAD_CRC entries
            continue;
          }

          norm[side][4] += 1.0L; // CRC-good entries in selected BCO ranges
          const ULong64_t frame_slot = bco / FRAME_BCO;
          observed_slots[side].insert(frame_slot);
          const int frame = static_cast<int>(frame_slot % N_FRAME_BINS);
          const double zabs = (static_cast<double>(frame) + 0.5) * Z_MAX_MM / N_FRAME_BINS;
          const double z_mm = side == 0 ? -zabs : zabs;

          for (int i = 0; i < 8; ++i)
          {
            const ULong64_t nsamples = t->dc_nsamples[i];
            if (nsamples == 0)
            {
              rejection[side][1] += 1.0L;
              continue;
            }

            const int channel = static_cast<int>(t->dc_channel) - (7 - i);
            const auto popt = get_pad(pad_lookup, static_cast<int>(t->dc_fee), channel);
            if (!popt || !popt->valid)
            {
              rejection[side][2] += 1.0L;
              continue;
            }
            const PadInfo& p = *popt;

            const double phi = wrap_phi(p.phi_base + (ebdc - 12) * M_PI / 6.0);
            const double raw = static_cast<double>(t->dc_current[i]);
            const double ped = raw - FIXED_PEDESTAL * static_cast<double>(nsamples);
            const double pedpos = std::max(0.0, ped);
            const bool strict_ok = ped >= 0.0 && ped <= STRICT_MAX_PEDSUB;

            // Exposure is stored independently so any variant can later be normalized.
            H.exposure_frame[side]->Fill(phi, p.r_mm, frame, static_cast<double>(nsamples));
            H.count_frame[side]->Fill(phi, p.r_mm, frame, 1.0);
            if (p.local_layer >= 0 && p.local_layer < 16 &&
                p.module >= 0 && p.module < N_MODULES &&
                p.local_pad >= 0 && p.local_pad < N_PADS[p.module])
            {
              H.native_exposure[side][sector][p.module]->Fill(
                  p.local_pad, p.local_layer, static_cast<double>(nsamples));
              H.native_count[side][sector][p.module]->Fill(p.local_pad, p.local_layer, 1.0);
            }

            fill_variant(H, RAW, side, sector, p, frame, z_mm, raw);
            fill_variant(H, PED60_SIGNED, side, sector, p, frame, z_mm, ped);
            fill_variant(H, PED60_POSITIVE, side, sector, p, frame, z_mm, pedpos);
            if (strict_ok)
              fill_variant(H, PED60_STRICT, side, sector, p, frame, z_mm, ped);

            norm[side][5] += 1.0L; // mapped sample blocks with nsamples > 0
            norm[side][6] += static_cast<long double>(nsamples);
            norm[side][7] += static_cast<long double>(raw);
            norm[side][8] += static_cast<long double>(ped);
            norm[side][9] += static_cast<long double>(pedpos);

            if (ped < 0.0) rejection[side][3] += 1.0L;
            if (ped > STRICT_MAX_PEDSUB) rejection[side][4] += 1.0L;
            if (strict_ok)
            {
              norm[side][10] += static_cast<long double>(ped);
              rejection[side][5] += 1.0L;
            }
          }
        }
      }
    }

    if (covered_frames != live_frame_slots)
    {
      std::cerr << "WARNING: indexed coverage for ebdc=" << ebdc
                << " stream=" << stream
                << " is " << covered_frames << " frames, expected "
                << live_frame_slots << std::endl;
    }
  }

  for (int side = 0; side < N_SIDES; ++side)
    norm[side][11] = static_cast<long double>(observed_slots[side].size());

  TFile fout(outputRoot, "RECREATE");
  if (fout.IsZombie())
  {
    std::cerr << "ERROR: cannot create " << outputRoot << std::endl;
    return;
  }

  for (int v = 0; v < N_VARIANTS; ++v)
  {
    for (int side = 0; side < N_SIDES; ++side)
    {
      H.frame[v][side]->Write();
      H.z[v][side]->Write();
    }
  }

  for (int side = 0; side < N_SIDES; ++side)
  {
    H.exposure_frame[side]->Write();
    H.count_frame[side]->Write();
    for (int sec = 0; sec < N_SECTORS; ++sec)
    {
      for (int mod = 0; mod < N_MODULES; ++mod)
      {
        for (int v = 0; v < N_VARIANTS; ++v)
          H.native[v][side][sec][mod]->Write();
        H.native_exposure[side][sec][mod]->Write();
        H.native_count[side][sec][mod]->Write();
      }
    }
  }

  // Compatibility aliases for the existing notebook. For the first pass use
  // the same cleaning as the existing T_DigitalCurrent.C: current - 60*nsamples,
  // accepting 0 <= pedestal-subtracted value <= 30000.
  const char* QA_NAME = "PHGarfieldRawHitsQA";
  for (int side = 0; side < N_SIDES; ++side)
  {
    std::unique_ptr<TH3D> h3(static_cast<TH3D*>(H.frame[PED60_STRICT][side]->Clone(
        Form("h_adc_phi_vs_radius_vs_tbin_side%d_%s", side, QA_NAME))));
    h3->SetTitle(Form("DC strict pedestal-subtracted input side %d;#phi;r [mm];DC frame mod 80", side));
    h3->Write();

    for (int sec = 0; sec < N_SECTORS; ++sec)
    {
      for (int mod = 0; mod < N_MODULES; ++mod)
      {
        std::unique_ptr<TH2D> h2(static_cast<TH2D*>(H.native[PED60_STRICT][side][sec][mod]->Clone(
            Form("h_adc_pad_vs_layer_side%d_sec%02d_R%d_%s", side, sec, mod+1, QA_NAME))));
        h2->Write();
      }
    }
  }

  TH1D hEvents(Form("h_nEvents_%s", QA_NAME),
               "DC equivalent TPC-volume exposure;quantity;count", 1, 0.5, 1.5);
  hEvents.GetXaxis()->SetBinLabel(1, "live_frame_slots/80");
  hEvents.SetBinContent(1, static_cast<double>(live_frame_slots) / 80.0);
  hEvents.Write();

  TH2D hNorm("h_dc_norm",
             "Additive DC normalization;TPC side;quantity",
             2, -0.5, 1.5, 16, 0.5, 16.5);
  hNorm.GetXaxis()->SetBinLabel(1, "South side0");
  hNorm.GetXaxis()->SetBinLabel(2, "North side1");
  const char* norm_labels[16] = {
      "N_WINDOWS", "LIVE_BCO_TICKS", "LIVE_FRAME_SLOTS", "LIVE_TPC_VOLUMES",
      "CRC_GOOD_ENTRIES", "MAPPED_SAMPLE_BLOCKS", "NSAMPLES_SUM", "RAW_CURRENT_SUM",
      "PED60_SIGNED_SUM", "PED60_POSITIVE_SUM", "PED60_STRICT_SUM", "OBSERVED_UNIQUE_FRAME_SLOTS",
      "GEM_R1_CURRENT_SUM", "GEM_R1_CURRENT2_SUM",
      "GEM_R1_CURRENT_X_LIVE_FRAMES", "GEM_R1_CURRENT_SAMPLES"};
  for (int m = 0; m < 16; ++m) hNorm.GetYaxis()->SetBinLabel(m+1, norm_labels[m]);
  for (int side = 0; side < N_SIDES; ++side)
    for (int m = 0; m < 16; ++m)
      hNorm.SetBinContent(side+1, m+1, static_cast<double>(norm[side][m]));
  hNorm.Write();

  // Convenient mergeable view. After hadd the bin contents are the
  // live-frame-weighted average R1 GEM load currents for the merged windows.
  TProfile pGem("p_dc_gem_r1_current",
                "R1 GEM load current, live-frame weighted;TPC side;load current",
                2, -0.5, 1.5);
  pGem.GetXaxis()->SetBinLabel(1, "South SR1");
  pGem.GetXaxis()->SetBinLabel(2, "North NR1");
  if (gem.available && live_frame_slots > 0)
  {
    pGem.Fill(0.0, gem.load_sr1, static_cast<double>(live_frame_slots));
    pGem.Fill(1.0, gem.load_nr1, static_cast<double>(live_frame_slots));
  }
  pGem.Write();

  TH2D hReject("h_dc_rejection_stats",
               "DC rejection/selection counters;TPC side;quantity",
               2, -0.5, 1.5, 6, 0.5, 6.5);
  hReject.GetXaxis()->SetBinLabel(1, "South side0");
  hReject.GetXaxis()->SetBinLabel(2, "North side1");
  const char* reject_labels[6] = {
      "BAD_CRC_ENTRIES", "NSAMPLES_ZERO", "BAD_MAPPING",
      "PED60_NEGATIVE", "PED60_GT30000", "PED60_STRICT_ACCEPTED"};
  for (int m = 0; m < 6; ++m) hReject.GetYaxis()->SetBinLabel(m+1, reject_labels[m]);
  for (int side = 0; side < N_SIDES; ++side)
    for (int m = 0; m < 6; ++m)
      hReject.SetBinContent(side+1, m+1, static_cast<double>(rejection[side][m]));
  hReject.Write();

  // One metadata entry per window. hadd concatenates these entries.
  TTree meta("dc_metadata", "DC synchronized-window metadata");
  int meta_run = run;
  int meta_window = window;
  int meta_nsub = static_cast<int>(win.size());
  ULong64_t meta_bco_start = bco_start;
  ULong64_t meta_bco_end = bco_end;
  ULong64_t meta_live_bco = live_bco_ticks;
  ULong64_t meta_live_frames = live_frame_slots;
  double meta_live_volumes = static_cast<double>(live_frame_slots) / 80.0;
  bool meta_conditions_available = gem.available;
  ULong64_t meta_conditions_bco = conditions_bco;
  ULong64_t meta_conditions_nearest_bco = gem.nearest_bco;
  ULong64_t meta_conditions_first_bco = gem.first_conditions_bco;
  ULong64_t meta_conditions_last_bco = gem.last_conditions_bco;
  float meta_load_sr1 = gem.load_sr1;
  float meta_load_nr1 = gem.load_nr1;
  float meta_average_sr1 = gem.average_sr1;
  float meta_average_nr1 = gem.average_nr1;
  std::string meta_cdb_global_tag = cdbGlobalTag;
  std::string meta_tpc_conditions_url = gem.cdb_url;
  meta.Branch("run", &meta_run);
  meta.Branch("window", &meta_window);
  meta.Branch("n_subintervals", &meta_nsub);
  meta.Branch("bco_start", &meta_bco_start);
  meta.Branch("bco_end", &meta_bco_end);
  meta.Branch("live_bco_ticks", &meta_live_bco);
  meta.Branch("live_frame_slots", &meta_live_frames);
  meta.Branch("live_tpc_volumes", &meta_live_volumes);
  meta.Branch("conditions_available", &meta_conditions_available);
  meta.Branch("conditions_bco", &meta_conditions_bco);
  meta.Branch("conditions_nearest_bco", &meta_conditions_nearest_bco);
  meta.Branch("conditions_first_bco", &meta_conditions_first_bco);
  meta.Branch("conditions_last_bco", &meta_conditions_last_bco);
  meta.Branch("gem_load_SR1", &meta_load_sr1);
  meta.Branch("gem_load_NR1", &meta_load_nr1);
  meta.Branch("gem_average_SR1", &meta_average_sr1);
  meta.Branch("gem_average_NR1", &meta_average_nr1);
  meta.Branch("cdb_global_tag", &meta_cdb_global_tag);
  meta.Branch("tpc_conditions_url", &meta_tpc_conditions_url);
  meta.Fill();
  meta.Write();

  TNamed("dc_pedestal_definition", "ped60 = current - 60*nsamples").Write();
  TNamed("dc_strict_definition", "ped60_strict accepts 0 <= current - 60*nsamples <= 30000").Write();
  TNamed("dc_frame_phase", "frame = floor(gtm_bco/9720) mod 80; validate absolute z phase before physics use").Write();
  TNamed("dc_gain_correction", "none; gain correction is intentionally deferred to the downstream notebook").Write();
  TNamed("dc_selected_notebook_variant", "ped60_strict").Write();
  TNamed("dc_gem_current_definition",
         "one R1 GEM load-current value per synchronized output file, evaluated at the live-exposure midpoint BCO; South=LoadSR1, North=LoadNR1; interpolation follows TpcConditionsReco").Write();

  fout.Close();
  std::cout << "\nWrote " << outputRoot << std::endl;
  std::cout << "Notebook compatibility variant: ped60_strict" << std::endl;
  std::cout << "Raw, signed pedestal-subtracted, positive-only, and strict maps are all retained." << std::endl;
}
