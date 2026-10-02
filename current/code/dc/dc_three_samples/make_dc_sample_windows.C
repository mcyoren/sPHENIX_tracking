#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dcsamples
{
  constexpr ULong64_t FRAME_BCO = 9720ULL;

  struct Interval
  {
    ULong64_t begin = 0;  // inclusive
    ULong64_t end = 0;    // exclusive
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

  std::vector<std::string> split(const std::string& line, char sep = ',')
  {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, sep)) out.push_back(item);
    return out;
  }

  std::vector<FileRecord> read_index(const std::string& filename, int run)
  {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open index: " + filename);

    std::vector<FileRecord> out;
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line))
    {
      if (line.empty()) continue;
      const auto f = split(line);
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

  std::vector<Interval> merge_intervals(std::vector<Interval> in)
  {
    if (in.empty()) return {};
    std::sort(in.begin(), in.end(), [](const Interval& a, const Interval& b)
    {
      if (a.begin != b.begin) return a.begin < b.begin;
      return a.end < b.end;
    });

    std::vector<Interval> out;
    out.push_back(in.front());
    for (size_t i = 1; i < in.size(); ++i)
    {
      if (in[i].begin <= out.back().end)
      {
        out.back().end = std::max(out.back().end, in[i].end);
      }
      else
      {
        out.push_back(in[i]);
      }
    }
    return out;
  }

  std::vector<Interval> intersect_sets(const std::vector<Interval>& a,
                                       const std::vector<Interval>& b)
  {
    std::vector<Interval> out;
    size_t i = 0;
    size_t j = 0;
    while (i < a.size() && j < b.size())
    {
      const ULong64_t lo = std::max(a[i].begin, b[j].begin);
      const ULong64_t hi = std::min(a[i].end, b[j].end);
      if (hi > lo) out.push_back({lo, hi});

      if (a[i].end < b[j].end) ++i;
      else ++j;
    }
    return merge_intervals(out);
  }

  ULong64_t interval_frames(const Interval& x)
  {
    if (x.end <= x.begin) return 0;
    return (x.end - x.begin) / FRAME_BCO;
  }

  // Find the common interval containing a target frame in the flattened
  // live-exposure coordinate. Return the interval index and local frame index.
  bool locate_flat_frame(const std::vector<Interval>& common,
                         ULong64_t flat_frame,
                         size_t& interval_index,
                         ULong64_t& local_frame)
  {
    ULong64_t cumulative = 0;
    for (size_t i = 0; i < common.size(); ++i)
    {
      const ULong64_t nf = interval_frames(common[i]);
      if (flat_frame < cumulative + nf)
      {
        interval_index = i;
        local_frame = flat_frame - cumulative;
        return true;
      }
      cumulative += nf;
    }
    return false;
  }
}

void make_dc_sample_windows(
    int run = 81558,
    ULong64_t sampleFrames = 2000,
    const char* indexCsv = "dc_index_run81558.csv",
    const char* outputIntervalsCsv = "dc_sample_windows_run81558.csv",
    const char* outputSummaryCsv = "dc_sample_windows_summary_run81558.csv",
    const char* outputJobs = "dc_sample_windows_run81558.jobs")
{
  using namespace dcsamples;

  if (sampleFrames == 0)
  {
    std::cerr << "ERROR: sampleFrames must be positive" << std::endl;
    return;
  }

  if (sampleFrames % 80 != 0)
  {
    std::cout << "WARNING: sampleFrames=" << sampleFrames
              << " is not a multiple of 80. For clean TPC-volume exposure "
              << "a multiple of 80 is preferable." << std::endl;
  }

  std::vector<FileRecord> files;
  try
  {
    files = read_index(indexCsv, run);
  }
  catch (const std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return;
  }

  std::array<std::vector<Interval>, 48> coverage;
  for (const auto& f : files)
  {
    if (f.ebdc < 0 || f.ebdc >= 24 || f.stream < 0 || f.stream >= 2) continue;
    const int key = 2 * f.ebdc + f.stream;

    // The last observed BCO represents a complete DC frame.
    coverage[key].push_back({f.first_bco, f.last_bco + FRAME_BCO});
  }

  for (int key = 0; key < 48; ++key)
  {
    coverage[key] = merge_intervals(coverage[key]);
    if (coverage[key].empty())
    {
      std::cerr << "ERROR: no coverage for ebdc=" << key / 2
                << " stream=" << key % 2 << std::endl;
      return;
    }
  }

  std::vector<Interval> common = coverage[0];
  for (int key = 1; key < 48; ++key)
  {
    common = intersect_sets(common, coverage[key]);
    if (common.empty())
    {
      std::cerr << "ERROR: common coverage becomes empty after ebdc="
                << key / 2 << " stream=" << key % 2 << std::endl;
      return;
    }
  }

  ULong64_t totalFrames = 0;
  for (const auto& x : common) totalFrames += interval_frames(x);

  if (totalFrames < sampleFrames)
  {
    std::cerr << "ERROR: total common exposure has only " << totalFrames
              << " frames, less than requested sampleFrames=" << sampleFrames
              << std::endl;
    return;
  }

  // Use 10%, 50%, and 90% of the flattened common live exposure.
  // This avoids startup/end boundaries while still sampling early/middle/late.
  const std::array<double, 3> fractions = {0.10, 0.50, 0.90};
  const std::array<const char*, 3> labels = {"early", "middle", "late"};

  struct Sample
  {
    int window = -1;
    std::string label;
    double fraction = 0.0;
    ULong64_t begin = 0;
    ULong64_t end = 0;
    ULong64_t frames = 0;
    size_t common_interval = 0;
  };

  std::vector<Sample> samples;
  samples.reserve(3);

  for (int w = 0; w < 3; ++w)
  {
    const ULong64_t targetFlat =
        static_cast<ULong64_t>(std::llround(fractions[w] * static_cast<double>(totalFrames - 1)));

    size_t ic = 0;
    ULong64_t localTarget = 0;
    if (!locate_flat_frame(common, targetFlat, ic, localTarget))
    {
      std::cerr << "ERROR: could not locate target sample " << labels[w] << std::endl;
      return;
    }

    const ULong64_t nf = interval_frames(common[ic]);
    if (nf < sampleFrames)
    {
      // Very unlikely for these data. Pick the largest common interval instead
      // rather than creating a sample across a gap.
      size_t best = ic;
      ULong64_t bestFrames = nf;
      for (size_t k = 0; k < common.size(); ++k)
      {
        const ULong64_t nfk = interval_frames(common[k]);
        if (nfk > bestFrames)
        {
          best = k;
          bestFrames = nfk;
        }
      }
      if (bestFrames < sampleFrames)
      {
        std::cerr << "ERROR: no single common interval is long enough for "
                  << sampleFrames << " frames" << std::endl;
        return;
      }
      ic = best;
      localTarget = bestFrames / 2;
    }

    const ULong64_t intervalFrames = interval_frames(common[ic]);
    ULong64_t startLocal = 0;
    if (localTarget > sampleFrames / 2)
      startLocal = localTarget - sampleFrames / 2;

    if (startLocal + sampleFrames > intervalFrames)
      startLocal = intervalFrames - sampleFrames;

    const ULong64_t begin = common[ic].begin + startLocal * FRAME_BCO;
    const ULong64_t end = begin + sampleFrames * FRAME_BCO;

    samples.push_back({w, labels[w], fractions[w], begin, end, sampleFrames, ic});
  }

  // Ensure the three samples are distinct/non-overlapping.
  for (size_t i = 1; i < samples.size(); ++i)
  {
    if (samples[i].begin < samples[i - 1].end)
    {
      std::cerr << "ERROR: requested sample width is too large; "
                << samples[i - 1].label << " and " << samples[i].label
                << " overlap." << std::endl;
      return;
    }
  }

  std::ofstream out(outputIntervalsCsv);
  out << "run,window,subinterval,bco_start,bco_end,live_frame_slots\n";
  for (const auto& s : samples)
  {
    out << run << ',' << s.window << ",0,"
        << s.begin << ',' << s.end << ',' << s.frames << '\n';
  }
  out.close();

  std::ofstream summary(outputSummaryCsv);
  summary << "run,window,label,fraction,bco_start,bco_end,n_subintervals,"
             "live_bco_ticks,live_frame_slots,live_tpc_volumes\n";
  for (const auto& s : samples)
  {
    summary << run << ',' << s.window << ',' << s.label << ','
            << s.fraction << ',' << s.begin << ',' << s.end << ",1,"
            << s.frames * FRAME_BCO << ',' << s.frames << ','
            << static_cast<double>(s.frames) / 80.0 << '\n';

    std::cout << "sample " << s.window << " (" << s.label << ", "
              << static_cast<int>(100.0 * s.fraction) << "%): "
              << s.frames << " frames = "
              << static_cast<double>(s.frames) / 80.0
              << " TPC volumes, BCO [" << s.begin << "," << s.end << ")"
              << std::endl;
  }
  summary.close();

  std::ofstream jobs(outputJobs);
  for (int w = 0; w < 3; ++w) jobs << w << '\n';
  jobs.close();

  std::cout << "\nCommon live exposure: " << totalFrames << " frames"
            << "\nEach sample: " << sampleFrames << " frames"
            << "\nTotal processed live exposure from 3 samples: "
            << 3 * sampleFrames << " frames ("
            << 100.0 * 3.0 * static_cast<double>(sampleFrames) /
                   static_cast<double>(totalFrames)
            << "% of common live exposure)"
            << "\n\nWrote:\n  " << outputIntervalsCsv
            << "\n  " << outputSummaryCsv
            << "\n  " << outputJobs << std::endl;
}
