#include <TString.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace dcwindows
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

  std::vector<std::string> split(const std::string& line, char sep=',')
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
    std::getline(in, line); // header
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
      // Half-open intervals that touch are continuous.
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
    size_t i = 0, j = 0;
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
}

void make_dc_windows(
    int run = 81558,
    int nWindows = 10,
    const char* indexCsv = "dc_index_run81558.csv",
    const char* outputIntervalsCsv = "dc_windows_run81558.csv",
    const char* outputSummaryCsv = "dc_windows_summary_run81558.csv",
    const char* outputJobs = "dc_windows_run81558.jobs")
{
  using namespace dcwindows;

  if (nWindows <= 0)
  {
    std::cerr << "ERROR: nWindows must be positive" << std::endl;
    return;
  }

  std::vector<FileRecord> files;
  try { files = read_index(indexCsv, run); }
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

    // Treat the last observed DC BCO as one complete 9720-BCO frame.
    coverage[key].push_back({f.first_bco, f.last_bco + FRAME_BCO});
  }

  for (int key = 0; key < 48; ++key)
  {
    coverage[key] = merge_intervals(coverage[key]);
    if (coverage[key].empty())
    {
      std::cerr << "ERROR: no coverage for ebdc=" << key/2
                << " stream=" << key%2 << std::endl;
      return;
    }
  }

  // Intersection of the union coverage of all 48 streams.
  std::vector<Interval> common = coverage[0];
  for (int key = 1; key < 48; ++key)
  {
    common = intersect_sets(common, coverage[key]);
    if (common.empty())
    {
      std::cerr << "ERROR: common coverage becomes empty after ebdc=" << key/2
                << " stream=" << key%2 << std::endl;
      return;
    }
  }

  ULong64_t totalFrames = 0;
  for (const auto& x : common)
  {
    const ULong64_t span = x.end - x.begin;
    if (span % FRAME_BCO != 0)
    {
      std::cerr << "WARNING: common interval is not an integer number of 9720-BCO frames: "
                << x.begin << " " << x.end << std::endl;
    }
    totalFrames += interval_frames(x);
  }

  if (totalFrames < static_cast<ULong64_t>(nWindows))
  {
    std::cerr << "ERROR: only " << totalFrames
              << " common frame slots for " << nWindows << " windows" << std::endl;
    return;
  }

  std::cout << "Common coverage has " << common.size() << " interval(s), "
            << totalFrames << " frame slots (" << totalFrames * FRAME_BCO
            << " BCO ticks)." << std::endl;

  // Desired number of live frame slots per window, distributing remainder first.
  std::vector<ULong64_t> target(nWindows, totalFrames / nWindows);
  for (int w = 0; w < totalFrames % static_cast<ULong64_t>(nWindows); ++w)
    target[w]++;

  struct Piece { int window; int sub; ULong64_t begin; ULong64_t end; ULong64_t frames; };
  std::vector<Piece> pieces;

  size_t ic = 0;
  ULong64_t pos = common[0].begin;

  for (int w = 0; w < nWindows; ++w)
  {
    ULong64_t need = target[w];
    int sub = 0;

    while (need > 0)
    {
      while (ic < common.size() && pos >= common[ic].end)
      {
        ++ic;
        if (ic < common.size()) pos = common[ic].begin;
      }
      if (ic >= common.size())
      {
        std::cerr << "ERROR: internal window partition ran out of common coverage" << std::endl;
        return;
      }

      const ULong64_t avail = (common[ic].end - pos) / FRAME_BCO;
      if (avail == 0)
      {
        pos = common[ic].end;
        continue;
      }

      const ULong64_t take = std::min(need, avail);
      const ULong64_t end = pos + take * FRAME_BCO;
      pieces.push_back({w, sub++, pos, end, take});
      pos = end;
      need -= take;
    }
  }

  std::ofstream out(outputIntervalsCsv);
  out << "run,window,subinterval,bco_start,bco_end,live_frame_slots\n";
  for (const auto& p : pieces)
  {
    out << run << ',' << p.window << ',' << p.sub << ','
        << p.begin << ',' << p.end << ',' << p.frames << '\n';
  }
  out.close();

  std::ofstream summary(outputSummaryCsv);
  summary << "run,window,bco_start,bco_end,n_subintervals,live_bco_ticks,live_frame_slots,live_tpc_volumes\n";

  for (int w = 0; w < nWindows; ++w)
  {
    ULong64_t begin = 0, end = 0, frames = 0;
    int nsub = 0;
    for (const auto& p : pieces)
    {
      if (p.window != w) continue;
      if (nsub == 0) begin = p.begin;
      end = p.end;
      frames += p.frames;
      ++nsub;
    }

    summary << run << ',' << w << ',' << begin << ',' << end << ','
            << nsub << ',' << frames * FRAME_BCO << ',' << frames << ','
            << static_cast<double>(frames) / 80.0 << '\n';

    std::cout << "window " << w
              << ": frames=" << frames
              << " subintervals=" << nsub
              << " BCO span " << begin << " -> " << end << std::endl;
  }
  summary.close();

  std::ofstream jobs(outputJobs);
  for (int w = 0; w < nWindows; ++w) jobs << w << '\n';
  jobs.close();

  std::cout << "\nWrote:\n  " << outputIntervalsCsv
            << "\n  " << outputSummaryCsv
            << "\n  " << outputJobs << std::endl;
}
