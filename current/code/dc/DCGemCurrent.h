#ifndef DCGEMCURRENT_H
#define DCGEMCURRENT_H

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Standalone helper that reproduces the R1 load-current calculation used by
// TpcConditionsReco, without requiring a Fun4All node tree / GL1 packet.
//
// The GEM channel lists and the interpolation below are copied from the
// TpcConditionsReco logic.  This is intentional: the DC preprocessing already
// knows the run number and a representative BCO for its synchronized window.
namespace dcgem
{
  const std::vector<std::string> SR1 = {
      "S_01_R1_G4_IMon",
      "S_04_R1_G4_IMon",
      "S_05_R1_G4_IMon",
      "S_10_R1_G4_IMon"};

  const std::vector<std::string> NR1 = {
      "N_03_R1_G4_IMon",
      "N_04_R1_G4_IMon",
      "N_06_R1_G4_IMon",
      "N_07_R1_G4_IMon",
      "N_08_R1_G4_IMon",
      "N_10_R1_G4_IMon",
      "N_11_R1_G4_IMon",
      "N_12_R1_G4_IMon"};

  struct Result
  {
    bool available = false;
    uint64_t target_bco = 0;              // 40-bit DC/GTM BCO
    uint64_t aligned_target_bco = 0;      // same time on absolute CDB BCO axis
    uint64_t nearest_bco = 0;             // absolute CDB BCO
    uint64_t first_conditions_bco = 0;
    uint64_t last_conditions_bco = 0;
    float load_sr1 = 0.0F;
    float load_nr1 = 0.0F;
    float average_sr1 = 0.0F;
    float average_nr1 = 0.0F;
    std::string cdb_url;
    std::string global_tag;
  };

  inline float median_current(CDBTTree& tree,
                              int channel,
                              const std::vector<std::string>& channels)
  {
    std::vector<float> currents;
    currents.reserve(channels.size());
    for (const auto& name : channels)
    {
      currents.push_back(tree.GetFloatValue(channel, name));
    }
    std::sort(currents.begin(), currents.end());
    const size_t n = currents.size();
    if (n == 0) return 0.0F;
    if (n % 2) return currents[n / 2];
    return 0.5F * (currents[n / 2 - 1] + currents[n / 2]);
  }

  inline float average_median_current(CDBTTree& tree,
                                      const std::map<uint64_t, int>& bco_to_channel,
                                      const std::vector<std::string>& channels)
  {
    if (bco_to_channel.empty()) return 0.0F;
    double sum = 0.0;
    for (const auto& [bco, channel] : bco_to_channel)
    {
      (void) bco;
      sum += median_current(tree, channel, channels);
    }
    return static_cast<float>(sum / static_cast<double>(bco_to_channel.size()));
  }

  // Same monotonic cubic Hermite interpolation currently used by
  // TpcConditionsReco::get_InterpolatedMedianCurrent.
  inline float interpolated_median_current(CDBTTree& tree,
                                           const std::map<uint64_t, int>& bco_to_channel,
                                           uint64_t bco,
                                           const std::vector<std::string>& channels)
  {
    if (bco_to_channel.empty()) return 0.0F;
    if (bco_to_channel.size() == 1)
      return median_current(tree, bco_to_channel.begin()->second, channels);

    const auto first = bco_to_channel.begin();
    const auto last = std::prev(bco_to_channel.end());
    if (bco <= first->first) return median_current(tree, first->second, channels);
    if (bco >= last->first) return median_current(tree, last->second, channels);

    const auto right = bco_to_channel.upper_bound(bco);
    const auto left = std::prev(right);

    const double x1 = static_cast<double>(left->first);
    const double x2 = static_cast<double>(right->first);
    const double y1 = median_current(tree, left->second, channels);
    const double y2 = median_current(tree, right->second, channels);
    const double h1 = x2 - x1;
    const double d1 = (y2 - y1) / h1;

    double m1 = d1;
    double m2 = d1;

    if (left != bco_to_channel.begin())
    {
      const auto left2 = std::prev(left);
      const double x0 = static_cast<double>(left2->first);
      const double y0 = median_current(tree, left2->second, channels);
      const double h0 = x1 - x0;
      const double d0 = (y1 - y0) / h0;
      if (d0 * d1 > 0.0)
      {
        const double w1 = 2.0 * h1 + h0;
        const double w2 = h1 + 2.0 * h0;
        m1 = (w1 + w2) / (w1 / d0 + w2 / d1);
      }
      else m1 = 0.0;
    }

    const auto right2 = std::next(right);
    if (right2 != bco_to_channel.end())
    {
      const double x3 = static_cast<double>(right2->first);
      const double y3 = median_current(tree, right2->second, channels);
      const double h2 = x3 - x2;
      const double d2 = (y3 - y2) / h2;
      if (d1 * d2 > 0.0)
      {
        const double w1 = 2.0 * h2 + h1;
        const double w2 = h2 + 2.0 * h1;
        m2 = (w1 + w2) / (w1 / d1 + w2 / d2);
      }
      else m2 = 0.0;
    }

    const double t = (static_cast<double>(bco) - x1) / h1;
    const double t2 = t * t;
    const double t3 = t2 * t;
    return static_cast<float>(
        (2.0 * t3 - 3.0 * t2 + 1.0) * y1 +
        (t3 - 2.0 * t2 + t) * h1 * m1 +
        (-2.0 * t3 + 3.0 * t2) * y2 +
        (t3 - t2) * h1 * m2);
  }

  inline Result load_r1_currents(int run,
                                 uint64_t target_bco,
                                 const std::string& global_tag = "newcdbtag")
  {
    Result out;
    out.target_bco = target_bco;
    out.global_tag = global_tag;

    recoConsts* rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", run);
    rc->set_IntFlag("RUNSEGMENT", 0);
    rc->set_StringFlag("CDB_GLOBALTAG", global_tag);
    rc->set_uint64Flag("TIMESTAMP", static_cast<uint64_t>(run));

    out.cdb_url = CDBInterface::instance()->getUrl("TPC_CONDITIONS");
    if (out.cdb_url.empty())
    {
      std::cerr << "WARNING: TPC_CONDITIONS is unavailable for run " << run
                << " with CDB tag " << global_tag << std::endl;
      return out;
    }

    CDBTTree tree(out.cdb_url);
    tree.LoadCalibrations();

    std::map<uint64_t, int> bco_to_channel;
    const auto& entries = tree.GetUInt64EntryMap();
    for (unsigned int channel = 0; channel < entries.size(); ++channel)
    {
      const uint64_t bco = tree.GetUInt64Value(channel, "bco");
      bco_to_channel[bco] = static_cast<int>(channel);
    }

    if (bco_to_channel.empty())
    {
      std::cerr << "WARNING: TPC_CONDITIONS has no BCO entries for run " << run << std::endl;
      return out;
    }

    out.first_conditions_bco = bco_to_channel.begin()->first;
    out.last_conditions_bco = std::prev(bco_to_channel.end())->first;

    // The TPC/DC GTM clock is a 40-bit counter.  TPC_CONDITIONS stores an
    // unwrapped/absolute BCO built from wall-clock time.  Lift the 40-bit DC
    // BCO onto the absolute CDB axis by choosing the rollover number nearest
    // the center of this run's conditions history.
    constexpr uint64_t GTM_RANGE = (1ULL << 40U);
    constexpr uint64_t GTM_MASK = GTM_RANGE - 1ULL;
    const uint64_t target40 = target_bco & GTM_MASK;
    const long double center = 0.5L *
        (static_cast<long double>(out.first_conditions_bco) +
         static_cast<long double>(out.last_conditions_bco));
    const long double k_real =
        (center - static_cast<long double>(target40)) /
        static_cast<long double>(GTM_RANGE);
    const int64_t k = static_cast<int64_t>(std::llround(k_real));
    const __int128 aligned128 = static_cast<__int128>(target40) +
                                static_cast<__int128>(k) *
                                static_cast<__int128>(GTM_RANGE);
    if (aligned128 < 0 || aligned128 > static_cast<__int128>(std::numeric_limits<uint64_t>::max()))
    {
      throw std::runtime_error("Failed to align 40-bit DC GTM BCO to TPC_CONDITIONS axis");
    }
    const uint64_t aligned_target = static_cast<uint64_t>(aligned128);
    out.aligned_target_bco = aligned_target;

    auto upper = bco_to_channel.upper_bound(aligned_target);
    auto selected = upper;
    if (upper == bco_to_channel.end())
    {
      selected = std::prev(bco_to_channel.end());
    }
    else if (upper == bco_to_channel.begin())
    {
      selected = upper;
    }
    else
    {
      const auto lower = std::prev(upper);
      const uint64_t dlo = aligned_target - lower->first;
      const uint64_t dhi = upper->first - aligned_target;
      selected = (dlo <= dhi) ? lower : upper;
    }
    out.nearest_bco = selected->first;

    out.load_sr1 = interpolated_median_current(tree, bco_to_channel, aligned_target, SR1);
    out.load_nr1 = interpolated_median_current(tree, bco_to_channel, aligned_target, NR1);
    out.average_sr1 = average_median_current(tree, bco_to_channel, SR1);
    out.average_nr1 = average_median_current(tree, bco_to_channel, NR1);
    out.available = true;

    std::cout << "TPC_CONDITIONS run " << run
              << " DC 40-bit BCO " << target40
              << " aligned absolute BCO " << aligned_target
              << " nearest stored BCO " << out.nearest_bco
              << " LoadSR1 " << out.load_sr1
              << " LoadNR1 " << out.load_nr1
              << std::endl;
    return out;
  }
}

#endif
