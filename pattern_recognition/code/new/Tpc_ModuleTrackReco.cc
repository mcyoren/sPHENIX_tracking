#include "Tpc_ModuleTrackReco.h"
#include "Tpc_FittingTools.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include "Tpc_ModuleTrack.h"
#include "Tpc_ModuleTrackContainer.h"
#include "Tpc_ModuleTrackContainerv1.h"
#include "Tpc_ModuleTrackv1.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <pthread.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <unordered_map>
#include <vector>

// ===================================================================
// Internal helpers (anonymous namespace)
// ===================================================================
namespace
{
  inline uint64_t make_blob_lookup_key(unsigned int layer,
                                       unsigned short pad,
                                       unsigned short tbin)
  {
    return (static_cast<uint64_t>(layer) << 32U) |
           (static_cast<uint64_t>(pad) << 16U) |
           static_cast<uint64_t>(tbin);
  }

  inline double sqr(const double value)
  {
    return value * value;
  }

  inline double safe_scale(const double value)
  {
    return std::max(std::fabs(value), 1.0e-6);
  }

  inline unsigned int module_first_layer(const InModuleThreadData* d)
  {
    return d->region * 16U + 7U;
  }

  inline unsigned int module_last_layer(const InModuleThreadData* d)
  {
    return module_first_layer(d) + 15U;
  }

  bool contains_blob(const std::vector<unsigned int>& chain, const unsigned int blob_index)
  {
    return std::find(chain.begin(), chain.end(), blob_index) != chain.end();
  }

  struct BlobAdcSort
  {
    const std::vector<InModuleThreadData::Blob>* blobs;
    explicit BlobAdcSort(const std::vector<InModuleThreadData::Blob>* b)
      : blobs(b)
    {
    }

    bool operator()(unsigned int a, unsigned int b) const
    {
      if ((*blobs)[a].adc != (*blobs)[b].adc)
      {
        return (*blobs)[a].adc > (*blobs)[b].adc;
      }
      return a < b;
    }
  };

  struct BlobInnerSort
  {
    const std::vector<InModuleThreadData::Blob>* blobs;
    explicit BlobInnerSort(const std::vector<InModuleThreadData::Blob>* b)
      : blobs(b)
    {
    }

    bool operator()(unsigned int a, unsigned int b) const
    {
      if ((*blobs)[a].layer != (*blobs)[b].layer)
      {
        return (*blobs)[a].layer < (*blobs)[b].layer;
      }
      if ((*blobs)[a].adc != (*blobs)[b].adc)
      {
        return (*blobs)[a].adc > (*blobs)[b].adc;
      }
      return a < b;
    }
  };

  // -------------------------------------------------------------------
  // Temporary weighted line fit. This intentionally keeps the same fit
  // definition as the original implementation, including ADC weighting.
  // -------------------------------------------------------------------
  bool fit_track_from_blobs(const std::vector<InModuleThreadData::Blob>& blobs,
                            const std::vector<unsigned int>& idx,
                            double weight_power,
                            double floor_frac,
                            InModuleThreadData::Track& trk)
  {
    if (idx.size() < 2)
    {
      return false;
    }

    double maxadc = 0.0;
    for (unsigned int i : idx)
    {
      maxadc = std::max(blobs[i].adc, maxadc);
    }

    std::vector<double> x;
    std::vector<double> pad;
    std::vector<double> tbin;
    std::vector<double> w;
    x.reserve(idx.size());
    pad.reserve(idx.size());
    tbin.reserve(idx.size());
    w.reserve(idx.size());

    unsigned int first_layer = std::numeric_limits<unsigned int>::max();
    unsigned int last_layer = 0;

    for (unsigned int i : idx)
    {
      const InModuleThreadData::Blob& bl = blobs[i];
      x.push_back(static_cast<double>(bl.layer));
      pad.push_back(bl.pad);
      tbin.push_back(bl.tbin);
      w.push_back(Tpc_FittingTools::adcWeight(bl.adc, maxadc, weight_power, floor_frac));

      first_layer = std::min(bl.layer, first_layer);
      last_layer = std::max(bl.layer, last_layer);
    }

    double mp = 0.0;
    double bp = 0.0;
    double cp = 0.0;
    double mt = 0.0;
    double bt = 0.0;
    double ct = 0.0;
    int ndp = 0;
    int ndt = 0;

    if (!Tpc_FittingTools::weightedLineFit(x, pad, w, mp, bp, cp, ndp))
    {
      return false;
    }
    if (!Tpc_FittingTools::weightedLineFit(x, tbin, w, mt, bt, ct, ndt))
    {
      return false;
    }

    trk.first_layer = first_layer;
    trk.last_layer = last_layer;
    trk.nblobs = static_cast<unsigned int>(idx.size());
    trk.pad_slope = mp;
    trk.pad_intercept = bp;
    trk.tbin_slope = mt;
    trk.tbin_intercept = bt;
    trk.blob_indices = idx;
    return true;
  }

  unsigned int count_raw_hits_from_blob_chain(const std::vector<InModuleThreadData::Blob>& blobs,
                                               const std::vector<unsigned int>& blob_idx)
  {
    unsigned int count = 0;
    for (unsigned int ib : blob_idx)
    {
      count += blobs[ib].raw_hit_count;
    }
    return count;
  }

  bool make_track_from_blob_chain(const std::vector<InModuleThreadData::Blob>& blobs,
                                  const std::vector<unsigned int>& blob_idx,
                                  double weight_power,
                                  double floor_frac,
                                  InModuleThreadData::Track& trk)
  {
    if (!fit_track_from_blobs(blobs, blob_idx, weight_power, floor_frac, trk))
    {
      return false;
    }

    trk.nblobs = static_cast<unsigned int>(blob_idx.size());
    trk.nrawhits = count_raw_hits_from_blob_chain(blobs, blob_idx);
    return trk.nrawhits > 0;
  }

  // -------------------------------------------------------------------
  // Track-piece connection
  // -------------------------------------------------------------------
  struct TrackStartSort
  {
    const std::vector<InModuleThreadData::Track>* tracks;
    explicit TrackStartSort(const std::vector<InModuleThreadData::Track>* t)
      : tracks(t)
    {
    }

    bool operator()(unsigned int a, unsigned int b) const
    {
      const InModuleThreadData::Track& ta = (*tracks)[a];
      const InModuleThreadData::Track& tb = (*tracks)[b];
      if (ta.first_layer != tb.first_layer)
      {
        return ta.first_layer < tb.first_layer;
      }
      if (ta.last_layer != tb.last_layer)
      {
        return ta.last_layer < tb.last_layer;
      }
      return ta.nrawhits > tb.nrawhits;
    }
  };

  void append_unique_blob_indices(std::vector<unsigned int>& dst,
                                  const std::vector<unsigned int>& src)
  {
    for (unsigned int i : src)
    {
      if (std::find(dst.begin(), dst.end(), i) == dst.end())
      {
        dst.push_back(i);
      }
    }
  }

  bool tracks_can_connect(const InModuleThreadData::Track& a,
                          const InModuleThreadData::Track& b,
                          unsigned int connect_max_layer_gap,
                          double connect_dp,
                          double connect_dt,
                          double connect_dpad_slope,
                          double connect_dtbin_slope,
                          double& score)
  {
    score = std::numeric_limits<double>::max();

    if (a.last_layer >= b.first_layer)
    {
      return false;
    }

    const unsigned int gap = b.first_layer - a.last_layer - 1;
    if (gap > connect_max_layer_gap)
    {
      return false;
    }

    const double lmatch = 0.5 * (static_cast<double>(a.last_layer) +
                                 static_cast<double>(b.first_layer));

    const double pad_a = a.pad_slope * lmatch + a.pad_intercept;
    const double pad_b = b.pad_slope * lmatch + b.pad_intercept;
    const double tbin_a = a.tbin_slope * lmatch + a.tbin_intercept;
    const double tbin_b = b.tbin_slope * lmatch + b.tbin_intercept;

    const double dp = std::fabs(pad_a - pad_b);
    const double dt = std::fabs(tbin_a - tbin_b);
    const double dmp = std::fabs(a.pad_slope - b.pad_slope);
    const double dmt = std::fabs(a.tbin_slope - b.tbin_slope);

    if (dp > connect_dp || dt > connect_dt ||
        dmp > connect_dpad_slope || dmt > connect_dtbin_slope)
    {
      return false;
    }

    score = sqr(dp / safe_scale(connect_dp)) +
            sqr(dt / safe_scale(connect_dt)) +
            sqr(dmp / safe_scale(connect_dpad_slope)) +
            sqr(dmt / safe_scale(connect_dtbin_slope)) +
            0.05 * static_cast<double>(gap);

    return true;
  }

  void connect_track_pieces_in_module(InModuleThreadData* d)
  {
    if (!d || d->tracks.size() < 2)
    {
      return;
    }

    std::vector<InModuleThreadData::Track> pieces = d->tracks;
    std::vector<InModuleThreadData::Track> output;
    std::vector<uint8_t> used(pieces.size(), 0);

    std::vector<unsigned int> order;
    order.reserve(pieces.size());
    for (unsigned int i = 0; i < pieces.size(); ++i)
    {
      order.push_back(i);
    }
    std::sort(order.begin(), order.end(), TrackStartSort(&pieces));

    for (unsigned int iseed : order)
    {
      if (used[iseed])
      {
        continue;
      }

      InModuleThreadData::Track current = pieces[iseed];
      used[iseed] = 1;

      bool merged_any = true;
      while (merged_any)
      {
        merged_any = false;
        int best_j = -1;
        double best_score = std::numeric_limits<double>::max();

        for (unsigned int j : order)
        {
          if (used[j])
          {
            continue;
          }

          double score = 0.0;
          if (!tracks_can_connect(current, pieces[j],
                                  d->connect_max_layer_gap,
                                  d->connect_dp,
                                  d->connect_dt,
                                  d->connect_dpad_slope,
                                  d->connect_dtbin_slope,
                                  score))
          {
            continue;
          }

          if (score < best_score)
          {
            best_score = score;
            best_j = static_cast<int>(j);
          }
        }

        if (best_j >= 0)
        {
          std::vector<unsigned int> trial = current.blob_indices;
          append_unique_blob_indices(trial, pieces[static_cast<unsigned int>(best_j)].blob_indices);
          std::sort(trial.begin(), trial.end(),
                    [d](unsigned int a, unsigned int b)
                    {
                      const auto& ba = d->blobs[a];
                      const auto& bb = d->blobs[b];
                      if (ba.layer != bb.layer)
                      {
                        return ba.layer < bb.layer;
                      }
                      return a < b;
                    });

          InModuleThreadData::Track refit = current;
          if (make_track_from_blob_chain(d->blobs, trial,
                                         d->weight_power,
                                         d->adc_weight_floor_frac,
                                         refit))
          {
            refit.pass = std::max(current.pass, pieces[static_cast<unsigned int>(best_j)].pass);
            refit.has_questionable = current.has_questionable || pieces[static_cast<unsigned int>(best_j)].has_questionable;
            refit.looper_candidate = current.looper_candidate || pieces[static_cast<unsigned int>(best_j)].looper_candidate;
            current = std::move(refit);
            used[static_cast<unsigned int>(best_j)] = 1;
            merged_any = true;
          }
        }
      }

      output.push_back(std::move(current));
    }

    for (unsigned int i = 0; i < output.size(); ++i)
    {
      output[i].track_id = i;
    }

    if (d->verbosity > 1)
    {
      std::cout << "Tpc_ModuleTrackReco connect pieces: region=" << d->region
                << " sector=" << d->sector << " side=" << d->side
                << " pieces=" << pieces.size()
                << " connected_tracks=" << output.size() << std::endl;
    }

    d->tracks.swap(output);
  }

  // -------------------------------------------------------------------
  // Conservative same-pad long-tail noise rejection. This retains exactly
  // the detector/electronics topology of the original cut, but replaces the
  // map-of-vectors temporary with one sorted index array.
  // -------------------------------------------------------------------
  struct RawHitLayerPadTimeSort
  {
    const std::vector<InModuleThreadData::RawHit>* raw_hits;
    explicit RawHitLayerPadTimeSort(const std::vector<InModuleThreadData::RawHit>* h)
      : raw_hits(h)
    {
    }

    bool operator()(unsigned int a, unsigned int b) const
    {
      const auto& ha = (*raw_hits)[a];
      const auto& hb = (*raw_hits)[b];
      if (ha.layer != hb.layer)
      {
        return ha.layer < hb.layer;
      }
      if (ha.pad != hb.pad)
      {
        return ha.pad < hb.pad;
      }
      if (ha.tbin != hb.tbin)
      {
        return ha.tbin < hb.tbin;
      }
      return ha.adc > hb.adc;
    }
  };

  void reject_long_pad_noise(InModuleThreadData* d)
  {
    if (!d || d->noise_max_consecutive_timebins <= 0 || d->raw_hits.empty())
    {
      return;
    }

    d->noise_keep_first_timebins = std::max(d->noise_keep_first_timebins, 0);

    std::vector<unsigned int> order(d->raw_hits.size());
    std::iota(order.begin(), order.end(), 0U);
    std::sort(order.begin(), order.end(), RawHitLayerPadTimeSort(&d->raw_hits));

    std::vector<uint8_t> remove(d->raw_hits.size(), 0);
    unsigned int nremoved = 0;

    unsigned int group_begin = 0;
    while (group_begin < order.size())
    {
      const auto& first = d->raw_hits[order[group_begin]];
      unsigned int group_end = group_begin + 1;
      while (group_end < order.size())
      {
        const auto& hit = d->raw_hits[order[group_end]];
        if (hit.layer != first.layer || hit.pad != first.pad)
        {
          break;
        }
        ++group_end;
      }

      if (group_end - group_begin > static_cast<unsigned int>(d->noise_max_consecutive_timebins))
      {
        unsigned int run_begin = group_begin;
        while (run_begin < group_end)
        {
          unsigned int run_end = run_begin;
          while (run_end + 1 < group_end)
          {
            const unsigned short t0 = d->raw_hits[order[run_end]].tbin;
            const unsigned short t1 = d->raw_hits[order[run_end + 1]].tbin;
            if (static_cast<int>(t1) != static_cast<int>(t0) + 1)
            {
              break;
            }
            ++run_end;
          }

          const unsigned int run_len = run_end - run_begin + 1;
          if (run_len > static_cast<unsigned int>(d->noise_max_consecutive_timebins))
          {
            const unsigned int keep_until =
                run_begin + static_cast<unsigned int>(d->noise_keep_first_timebins);

            for (unsigned int ir = run_begin; ir <= run_end; ++ir)
            {
              if (ir < keep_until || ir == run_begin)
              {
                continue;
              }

              const unsigned int cur_idx = order[ir];
              const unsigned int prev_idx = order[ir - 1];
              const int cur_adc = static_cast<int>(d->raw_hits[cur_idx].adc);
              const int prev_adc = static_cast<int>(d->raw_hits[prev_idx].adc);

              if (cur_adc <= prev_adc + d->noise_adc_tolerance)
              {
                if (!remove[cur_idx])
                {
                  remove[cur_idx] = 1;
                  ++nremoved;
                }
              }
            }
          }

          run_begin = run_end + 1;
        }
      }

      group_begin = group_end;
    }

    if (nremoved == 0)
    {
      return;
    }

    std::vector<InModuleThreadData::RawHit> kept;
    kept.reserve(d->raw_hits.size() - nremoved);
    for (unsigned int i = 0; i < d->raw_hits.size(); ++i)
    {
      if (!remove[i])
      {
        kept.push_back(d->raw_hits[i]);
      }
    }
    d->raw_hits.swap(kept);

    if (d->verbosity > 1)
    {
      std::cout << "Tpc_ModuleTrackReco noise rejection: region=" << d->region
                << " sector=" << d->sector << " side=" << d->side
                << " removed " << nremoved << " long same-pad tail hits"
                << std::endl;
    }
  }

  void collect_raw_hits(InModuleThreadData* d)
  {
    d->raw_hits.clear();

    for (unsigned int ihs = 0; ihs < d->layer_hitsets.size(); ++ihs)
    {
      TrkrHitSet* hitset = d->layer_hitsets[ihs].hitset;
      if (!hitset)
      {
        continue;
      }

      TrkrHitSet::ConstRange range = hitset->getHits();
      for (TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
      {
        const TrkrDefs::hitkey hitkey = hitr->first;
        const TrkrHit* hit = hitr->second;
        if (!hit)
        {
          continue;
        }

        const unsigned short pad = TpcDefs::getPad(hitkey);
        const unsigned short tbin = TpcDefs::getTBin(hitkey);
        const unsigned short rawAdc = hit->getAdc();
        const double fadc = static_cast<double>(rawAdc) - d->pedestal;
        if (fadc <= 0.0)
        {
          continue;
        }

        InModuleThreadData::RawHit rh;
        rh.layer = d->layer_hitsets[ihs].layer;
        rh.hitsetkey = d->layer_hitsets[ihs].hitsetkey;
        rh.hitkey = hitkey;
        rh.pad = pad;
        rh.tbin = tbin;
        rh.adc = static_cast<unsigned short>(fadc);
        d->raw_hits.push_back(rh);
      }
    }
  }

  void build_blob_layer_offsets(InModuleThreadData* d)
  {
    d->blob_layer_offsets.fill(static_cast<unsigned int>(d->blobs.size()));
    const unsigned int first_layer = module_first_layer(d);

    unsigned int iblob = 0;
    for (unsigned int ilocal = 0; ilocal < 16; ++ilocal)
    {
      const unsigned int layer = first_layer + ilocal;
      while (iblob < d->blobs.size() && d->blobs[iblob].layer < layer)
      {
        ++iblob;
      }
      d->blob_layer_offsets[ilocal] = iblob;
      while (iblob < d->blobs.size() && d->blobs[iblob].layer == layer)
      {
        ++iblob;
      }
      d->blob_layer_offsets[ilocal + 1] = iblob;
    }
  }

  void build_blobs(InModuleThreadData* d)
  {
    d->blobs.clear();
    d->blob_raw_hit_indices.clear();
    d->blob_layer_offsets.fill(0U);

    const unsigned int n = static_cast<unsigned int>(d->raw_hits.size());
    if (n == 0)
    {
      return;
    }

    std::unordered_map<uint64_t, unsigned int> hit_lookup;
    hit_lookup.reserve(n);

    for (unsigned int i = 0; i < n; ++i)
    {
      const InModuleThreadData::RawHit& hit = d->raw_hits[i];
      const uint64_t key = make_blob_lookup_key(hit.layer, hit.pad, hit.tbin);
      const auto inserted = hit_lookup.emplace(key, i);
      if (!inserted.second && d->verbosity > 2)
      {
        std::cout << "Tpc_ModuleTrackReco::build_blobs - duplicate "
                  << "(layer,pad,tbin)=(" << hit.layer << ","
                  << hit.pad << "," << hit.tbin << ")" << std::endl;
      }
    }

    std::vector<uint8_t> used(n, 0);
    d->blob_raw_hit_indices.reserve(n);

    for (unsigned int i = 0; i < n; ++i)
    {
      if (used[i])
      {
        continue;
      }

      used[i] = 1;
      std::deque<unsigned int> q;
      q.push_back(i);

      double sw = 0.0;
      double sp = 0.0;
      double st = 0.0;
      unsigned int nh = 0;
      const unsigned int layer = d->raw_hits[i].layer;

      InModuleThreadData::Blob bl;
      bl.raw_hit_begin = static_cast<unsigned int>(d->blob_raw_hit_indices.size());

      while (!q.empty())
      {
        const unsigned int a = q.front();
        q.pop_front();

        const InModuleThreadData::RawHit& ha = d->raw_hits[a];
        d->blob_raw_hit_indices.push_back(a);

        const double wa = static_cast<double>(ha.adc);
        sw += wa;
        sp += wa * static_cast<double>(ha.pad);
        st += wa * static_cast<double>(ha.tbin);
        ++nh;

        for (int dp = -d->blob_dp; dp <= d->blob_dp; ++dp)
        {
          const int neighbor_pad = static_cast<int>(ha.pad) + dp;
          if (neighbor_pad < 0 ||
              neighbor_pad > static_cast<int>(std::numeric_limits<unsigned short>::max()))
          {
            continue;
          }

          for (int dt = -d->blob_dt; dt <= d->blob_dt; ++dt)
          {
            const int neighbor_tbin = static_cast<int>(ha.tbin) + dt;
            if (neighbor_tbin < 0 ||
                neighbor_tbin > static_cast<int>(std::numeric_limits<unsigned short>::max()))
            {
              continue;
            }

            const uint64_t key = make_blob_lookup_key(
                layer,
                static_cast<unsigned short>(neighbor_pad),
                static_cast<unsigned short>(neighbor_tbin));
            const auto it = hit_lookup.find(key);
            if (it == hit_lookup.end())
            {
              continue;
            }

            const unsigned int j = it->second;
            if (used[j])
            {
              continue;
            }

            used[j] = 1;
            q.push_back(j);
          }
        }
      }

      if (sw <= 0.0)
      {
        d->blob_raw_hit_indices.resize(bl.raw_hit_begin);
        continue;
      }

      bl.layer = layer;
      bl.pad = sp / sw;
      bl.tbin = st / sw;
      bl.adc = sw;
      bl.nhits = nh;
      bl.raw_hit_count = static_cast<unsigned int>(d->blob_raw_hit_indices.size()) - bl.raw_hit_begin;
      bl.ownership = InModuleThreadData::BLOB_FREE;
      bl.owner_tracklet = -1;
      bl.owner_cost = std::numeric_limits<float>::infinity();
      d->blobs.push_back(bl);
    }

    // collect_raw_hits() visits module hitsets in increasing layer order and
    // the noise compaction preserves that order. Flood filling never crosses a
    // layer, so blobs are already layer-contiguous here. Do not reorder Blob
    // objects after filling blob_raw_hit_indices: raw_hit_begin/raw_hit_count
    // deliberately refer to that flat insertion order.
#ifndef NDEBUG
    for (unsigned int i = 1; i < d->blobs.size(); ++i)
    {
      if (d->blobs[i].layer < d->blobs[i - 1].layer)
      {
        std::cerr << "Tpc_ModuleTrackReco::build_blobs - blobs are not layer ordered"
                  << std::endl;
        break;
      }
    }
#endif
    build_blob_layer_offsets(d);
  }

  std::pair<unsigned int, unsigned int> blob_range_for_layer(const InModuleThreadData* d,
                                                             unsigned int layer)
  {
    const unsigned int first = module_first_layer(d);
    if (layer < first || layer > first + 15U)
    {
      return {0U, 0U};
    }
    const unsigned int local = layer - first;
    return {d->blob_layer_offsets[local], d->blob_layer_offsets[local + 1]};
  }

  // -------------------------------------------------------------------
  // Local first/second derivative information. Slopes are always defined
  // with increasing detector layer, even when the chain is being grown inward.
  // -------------------------------------------------------------------
  struct DerivativeState
  {
    bool first_valid{false};
    bool second_valid{false};

    double edge_pad_slope{0.0};
    double edge_tbin_slope{0.0};
    double edge_mid_layer{0.0};

    double mean_pad_slope{0.0};
    double mean_tbin_slope{0.0};
    double mean_mid_layer{0.0};
    double mean_pad_second{0.0};
    double mean_tbin_second{0.0};

    double local_pred_pad{0.0};
    double local_pred_tbin{0.0};
  };

  DerivativeState compute_derivative_state(const InModuleThreadData* d,
                                           const std::vector<unsigned int>& chain,
                                           int direction,
                                           unsigned int target_layer)
  {
    DerivativeState state;
    if (chain.size() < 2)
    {
      return state;
    }

    struct Segment
    {
      double midpoint{0.0};
      double pad_slope{0.0};
      double tbin_slope{0.0};
    };

    std::array<Segment, 16> segments{};
    unsigned int nseg = 0;
    for (unsigned int i = 1; i < chain.size() && nseg < segments.size(); ++i)
    {
      const auto& a = d->blobs[chain[i - 1]];
      const auto& b = d->blobs[chain[i]];
      const int dl_i = static_cast<int>(b.layer) - static_cast<int>(a.layer);
      if (dl_i <= 0)
      {
        continue;
      }
      const double dl = static_cast<double>(dl_i);
      segments[nseg].midpoint = 0.5 * (static_cast<double>(a.layer) + static_cast<double>(b.layer));
      segments[nseg].pad_slope = (b.pad - a.pad) / dl;
      segments[nseg].tbin_slope = (b.tbin - a.tbin) / dl;
      ++nseg;
    }

    if (nseg == 0)
    {
      return state;
    }

    const unsigned int history = std::max(1U, d->derivative_history);
    const unsigned int begin = direction > 0 ? (nseg > history ? nseg - history : 0U) : 0U;
    const unsigned int end = direction > 0 ? nseg : std::min(nseg, history);

    double sum_pad = 0.0;
    double sum_tbin = 0.0;
    double sum_mid = 0.0;
    for (unsigned int i = begin; i < end; ++i)
    {
      sum_pad += segments[i].pad_slope;
      sum_tbin += segments[i].tbin_slope;
      sum_mid += segments[i].midpoint;
    }
    const double inv = 1.0 / static_cast<double>(end - begin);
    state.mean_pad_slope = sum_pad * inv;
    state.mean_tbin_slope = sum_tbin * inv;
    state.mean_mid_layer = sum_mid * inv;
    state.first_valid = true;

    const unsigned int edge_index = direction > 0 ? end - 1U : begin;
    state.edge_pad_slope = segments[edge_index].pad_slope;
    state.edge_tbin_slope = segments[edge_index].tbin_slope;
    state.edge_mid_layer = segments[edge_index].midpoint;

    if (end - begin >= 2)
    {
      double sum_pad_second = 0.0;
      double sum_tbin_second = 0.0;
      unsigned int nsecond = 0;
      for (unsigned int i = begin + 1; i < end; ++i)
      {
        const double dmid = segments[i].midpoint - segments[i - 1].midpoint;
        if (std::fabs(dmid) < 1.0e-9)
        {
          continue;
        }
        sum_pad_second += (segments[i].pad_slope - segments[i - 1].pad_slope) / dmid;
        sum_tbin_second += (segments[i].tbin_slope - segments[i - 1].tbin_slope) / dmid;
        ++nsecond;
      }
      if (nsecond > 0)
      {
        state.mean_pad_second = sum_pad_second / static_cast<double>(nsecond);
        state.mean_tbin_second = sum_tbin_second / static_cast<double>(nsecond);
        state.second_valid = true;
      }
    }

    const auto& edge = direction > 0 ? d->blobs[chain.back()] : d->blobs[chain.front()];
    const double dl = static_cast<double>(static_cast<int>(target_layer) - static_cast<int>(edge.layer));
    state.local_pred_pad = edge.pad + state.edge_pad_slope * dl;
    state.local_pred_tbin = edge.tbin + state.edge_tbin_slope * dl;
    if (state.second_valid)
    {
      state.local_pred_pad += 0.5 * state.mean_pad_second * dl * dl;
      state.local_pred_tbin += 0.5 * state.mean_tbin_second * dl * dl;
    }

    return state;
  }

  struct CandidateEval
  {
    unsigned int blob_index{0};
    double straight_dp{0.0};
    double straight_dt{0.0};
    double local_dp{0.0};
    double local_dt{0.0};
    double first_pad_residual{0.0};
    double first_tbin_residual{0.0};
    double second_pad_residual{0.0};
    double second_tbin_residual{0.0};
    double cost{std::numeric_limits<double>::max()};
    double rank_score{std::numeric_limits<double>::max()};
    bool local_valid{false};
    bool first_valid{false};
    bool second_valid{false};
    bool accepted{false};
    int ownership_relation{0};  // 0 normal, 1 transfer, 2 share, 3 worse
  };

  void compute_candidate_derivative_residuals(const InModuleThreadData* d,
                                              const std::vector<unsigned int>& chain,
                                              int direction,
                                              const DerivativeState& deriv,
                                              CandidateEval& eval)
  {
    if (!deriv.first_valid || chain.empty())
    {
      return;
    }

    const auto& candidate = d->blobs[eval.blob_index];
    const auto& edge = direction > 0 ? d->blobs[chain.back()] : d->blobs[chain.front()];

    const InModuleThreadData::Blob* low = &edge;
    const InModuleThreadData::Blob* high = &candidate;
    if (candidate.layer < edge.layer)
    {
      low = &candidate;
      high = &edge;
    }

    const double dl = static_cast<double>(high->layer - low->layer);
    if (dl <= 0.0)
    {
      return;
    }

    const double candidate_mid = 0.5 * (static_cast<double>(low->layer) + static_cast<double>(high->layer));
    const double candidate_pad_slope = (high->pad - low->pad) / dl;
    const double candidate_tbin_slope = (high->tbin - low->tbin) / dl;

    // First derivative is compared to the recent average, propagated to the
    // candidate midpoint with the average second derivative when available.
    double expected_pad_slope = deriv.mean_pad_slope;
    double expected_tbin_slope = deriv.mean_tbin_slope;
    if (deriv.second_valid)
    {
      expected_pad_slope += deriv.mean_pad_second * (candidate_mid - deriv.mean_mid_layer);
      expected_tbin_slope += deriv.mean_tbin_second * (candidate_mid - deriv.mean_mid_layer);
    }

    eval.first_pad_residual = candidate_pad_slope - expected_pad_slope;
    eval.first_tbin_residual = candidate_tbin_slope - expected_tbin_slope;
    eval.first_valid = true;

    const double dmid = candidate_mid - deriv.edge_mid_layer;
    if (deriv.second_valid && std::fabs(dmid) > 1.0e-9)
    {
      const double candidate_pad_second = (candidate_pad_slope - deriv.edge_pad_slope) / dmid;
      const double candidate_tbin_second = (candidate_tbin_slope - deriv.edge_tbin_slope) / dmid;
      eval.second_pad_residual = candidate_pad_second - deriv.mean_pad_second;
      eval.second_tbin_residual = candidate_tbin_second - deriv.mean_tbin_second;
      eval.second_valid = true;
    }
  }

  bool derivative_cuts_pass(const InModuleThreadData* d,
                            unsigned int pass,
                            const CandidateEval& eval)
  {
    const double max_pad_first = pass == 1 ? d->pass1_max_dpad_slope_residual : d->pass2_max_dpad_slope_residual;
    const double max_tbin_first = pass == 1 ? d->pass1_max_dtbin_slope_residual : d->pass2_max_dtbin_slope_residual;
    const double max_pad_second = pass == 1 ? d->pass1_max_dpad_second_residual : d->pass2_max_dpad_second_residual;
    const double max_tbin_second = pass == 1 ? d->pass1_max_dtbin_second_residual : d->pass2_max_dtbin_second_residual;

    if (eval.first_valid)
    {
      if (max_pad_first >= 0.0 && std::fabs(eval.first_pad_residual) > max_pad_first)
      {
        return false;
      }
      if (max_tbin_first >= 0.0 && std::fabs(eval.first_tbin_residual) > max_tbin_first)
      {
        return false;
      }
    }

    if (eval.second_valid)
    {
      if (max_pad_second >= 0.0 && std::fabs(eval.second_pad_residual) > max_pad_second)
      {
        return false;
      }
      if (max_tbin_second >= 0.0 && std::fabs(eval.second_tbin_residual) > max_tbin_second)
      {
        return false;
      }
    }

    return true;
  }

  double association_cost(const InModuleThreadData* d,
                          unsigned int pass,
                          const CandidateEval& eval)
  {
    double position_cost = sqr(eval.straight_dp / safe_scale(d->score_dp_scale)) +
                           sqr(eval.straight_dt / safe_scale(d->score_dt_scale));

    if (pass == 2 && eval.local_valid)
    {
      const double local_cost = sqr(eval.local_dp / safe_scale(d->score_dp_scale)) +
                                sqr(eval.local_dt / safe_scale(d->score_dt_scale));
      position_cost = std::min(position_cost, local_cost);
    }

    double cost = position_cost;
    if (eval.first_valid)
    {
      cost += sqr(eval.first_pad_residual / safe_scale(d->score_dpad_slope_scale));
      cost += sqr(eval.first_tbin_residual / safe_scale(d->score_dtbin_slope_scale));
    }
    if (eval.second_valid)
    {
      cost += sqr(eval.second_pad_residual / safe_scale(d->score_dpad_second_scale));
      cost += sqr(eval.second_tbin_residual / safe_scale(d->score_dtbin_second_scale));
    }
    return cost;
  }

  int ownership_relation(const InModuleThreadData* d,
                         unsigned int track_index,
                         const InModuleThreadData::Blob& bl,
                         double new_cost)
  {
    if (bl.ownership == InModuleThreadData::BLOB_FREE)
    {
      return 0;
    }
    if (bl.ownership == InModuleThreadData::BLOB_HARD)
    {
      return bl.owner_tracklet == static_cast<int32_t>(track_index) ? 0 : 3;
    }
    if (!d->enable_questionable_reassignment)
    {
      return bl.owner_tracklet == static_cast<int32_t>(track_index) ? 0 : 3;
    }
    if (bl.ownership == InModuleThreadData::BLOB_SHARED)
    {
      return 2;
    }
    if (bl.owner_tracklet == static_cast<int32_t>(track_index))
    {
      return 0;
    }

    const double old_cost = static_cast<double>(bl.owner_cost);
    if (!std::isfinite(old_cost))
    {
      return 1;
    }
    // A nearly perfect existing match must not be stolen merely because its
    // stored cost is zero. Only another nearly perfect match is treated as
    // comparable/shared in that limit.
    if (old_cost <= 1.0e-9)
    {
      return new_cost <= 1.0e-9 ? 2 : 3;
    }
    if (new_cost < d->transfer_cost_ratio * old_cost)
    {
      return 1;
    }
    const double ratio = std::max(new_cost, old_cost) / std::max(std::min(new_cost, old_cost), 1.0e-9);
    if (ratio <= d->share_cost_ratio)
    {
      return 2;
    }
    return 3;
  }

  struct LayerSearchResult
  {
    unsigned int broad_count{0};
    std::vector<CandidateEval> candidates;
    int best_index{-1};
    double best_cost{std::numeric_limits<double>::max()};
    double second_cost{std::numeric_limits<double>::max()};
  };

  LayerSearchResult search_layer(const InModuleThreadData* d,
                                 unsigned int track_index,
                                 const std::vector<unsigned int>& chain,
                                 unsigned int target_layer,
                                 int direction,
                                 unsigned int pass)
  {
    LayerSearchResult result;
    if (chain.empty())
    {
      return result;
    }

    const auto& edge = direction > 0 ? d->blobs[chain.back()] : d->blobs[chain.front()];

    double straight_pred_pad = edge.pad;
    double straight_pred_tbin = edge.tbin;
    if (chain.size() >= 2)
    {
      InModuleThreadData::Track tmp;
      if (fit_track_from_blobs(d->blobs, chain,
                               d->weight_power,
                               d->adc_weight_floor_frac,
                               tmp))
      {
        straight_pred_pad = tmp.pad_slope * static_cast<double>(target_layer) + tmp.pad_intercept;
        straight_pred_tbin = tmp.tbin_slope * static_cast<double>(target_layer) + tmp.tbin_intercept;
      }
    }

    const DerivativeState deriv = compute_derivative_state(d, chain, direction, target_layer);
    const auto range = blob_range_for_layer(d, target_layer);
    result.candidates.reserve(range.second - range.first);

    for (unsigned int ib = range.first; ib < range.second; ++ib)
    {
      if (contains_blob(chain, ib))
      {
        continue;
      }

      const auto& bl = d->blobs[ib];
      if (bl.ownership == InModuleThreadData::BLOB_HARD &&
          bl.owner_tracklet != static_cast<int32_t>(track_index))
      {
        continue;
      }
      if (!d->enable_questionable_reassignment &&
          bl.ownership != InModuleThreadData::BLOB_FREE &&
          bl.owner_tracklet != static_cast<int32_t>(track_index))
      {
        continue;
      }

      CandidateEval eval;
      eval.blob_index = ib;
      eval.straight_dp = bl.pad - straight_pred_pad;
      eval.straight_dt = bl.tbin - straight_pred_tbin;
      eval.local_valid = deriv.first_valid;
      if (eval.local_valid)
      {
        eval.local_dp = bl.pad - deriv.local_pred_pad;
        eval.local_dt = bl.tbin - deriv.local_pred_tbin;
      }

      bool in_broad = false;
      if (pass == 1)
      {
        in_broad = std::fabs(eval.straight_dp) <= d->pass1_broad_dp &&
                   std::fabs(eval.straight_dt) <= d->pass1_broad_dt;
      }
      else
      {
        const bool straight_broad = std::fabs(eval.straight_dp) <= d->pass2_broad_dp &&
                                    std::fabs(eval.straight_dt) <= d->pass2_broad_dt;
        const bool local_broad = eval.local_valid &&
                                 std::fabs(eval.local_dp) <= d->pass2_broad_dp &&
                                 std::fabs(eval.local_dt) <= d->pass2_broad_dt;
        in_broad = straight_broad || local_broad;
      }

      if (!in_broad)
      {
        continue;
      }

      ++result.broad_count;
      compute_candidate_derivative_residuals(d, chain, direction, deriv, eval);

      bool position_ok = false;
      if (pass == 1)
      {
        position_ok = std::fabs(eval.straight_dp) <= d->pass1_tight_dp &&
                      std::fabs(eval.straight_dt) <= d->pass1_tight_dt;
      }
      else
      {
        const bool straight_ok = std::fabs(eval.straight_dp) <= d->pass2_straight_dp &&
                                 std::fabs(eval.straight_dt) <= d->pass2_straight_dt;
        const bool local_ok = eval.local_valid &&
                              std::fabs(eval.local_dp) <= d->pass2_local_dp &&
                              std::fabs(eval.local_dt) <= d->pass2_local_dt;
        const bool combined_ok = eval.local_valid &&
                                 std::fabs(eval.straight_dp) <= d->pass2_combined_dp &&
                                 std::fabs(eval.straight_dt) <= d->pass2_combined_dt &&
                                 std::fabs(eval.local_dp) <= d->pass2_combined_dp &&
                                 std::fabs(eval.local_dt) <= d->pass2_combined_dt;
        position_ok = straight_ok || local_ok || combined_ok;
      }

      eval.accepted = position_ok && derivative_cuts_pass(d, pass, eval);
      eval.cost = association_cost(d, pass, eval);
      eval.rank_score = eval.cost - 0.01 * std::log(bl.adc + 1.0);
      eval.ownership_relation = ownership_relation(d, track_index, bl, eval.cost);
      if (eval.ownership_relation == 3)
      {
        eval.accepted = false;
      }

      result.candidates.push_back(eval);
    }

    for (unsigned int i = 0; i < result.candidates.size(); ++i)
    {
      const auto& eval = result.candidates[i];
      if (!eval.accepted)
      {
        continue;
      }
      if (result.best_index < 0 ||
          eval.rank_score < result.candidates[static_cast<unsigned int>(result.best_index)].rank_score)
      {
        if (result.best_index >= 0)
        {
          result.second_cost = result.best_cost;
        }
        result.best_index = static_cast<int>(i);
        result.best_cost = eval.cost;
      }
      else if (eval.cost < result.second_cost)
      {
        result.second_cost = eval.cost;
      }
    }

    return result;
  }

  void record_layer_search_qa(InModuleThreadData* d,
                              const LayerSearchResult& result,
                              unsigned int pass,
                              unsigned int layer_step)
  {
    if (!d->do_pattern_qa)
    {
      return;
    }

    for (unsigned int i = 0; i < result.candidates.size(); ++i)
    {
      const auto& eval = result.candidates[i];
      InModuleThreadData::PatternQAEntry entry;
      entry.pass = static_cast<uint8_t>(pass);
      entry.selected = result.best_index == static_cast<int>(i) ? 1U : 0U;
      if (entry.selected)
      {
        const auto& selected = result.candidates[i];
        if (selected.ownership_relation == 2)
        {
          entry.association_state = 3U;  // shared
        }
        else if (selected.ownership_relation == 1)
        {
          entry.association_state = 4U;  // ownership transfer
        }
        else
        {
          entry.association_state = result.broad_count > 1U ? 2U : 1U;
        }
      }
      entry.layer_step = static_cast<uint8_t>(layer_step);
      entry.candidate_count = static_cast<uint8_t>(std::min(result.broad_count, 255U));
      entry.local_valid = eval.local_valid ? 1U : 0U;
      entry.first_valid = eval.first_valid ? 1U : 0U;
      entry.second_valid = eval.second_valid ? 1U : 0U;
      entry.straight_dp = static_cast<float>(eval.straight_dp);
      entry.straight_dt = static_cast<float>(eval.straight_dt);
      entry.local_dp = static_cast<float>(eval.local_dp);
      entry.local_dt = static_cast<float>(eval.local_dt);
      entry.first_pad_residual = static_cast<float>(eval.first_pad_residual);
      entry.first_tbin_residual = static_cast<float>(eval.first_tbin_residual);
      entry.second_pad_residual = static_cast<float>(eval.second_pad_residual);
      entry.second_tbin_residual = static_cast<float>(eval.second_tbin_residual);
      entry.cost = static_cast<float>(eval.cost);
      entry.second_cost = static_cast<float>(result.second_cost);
      d->pattern_qa.push_back(entry);
    }
  }

  void remove_questionable_records_from(InModuleThreadData* d,
                                        unsigned int track_index,
                                        unsigned int chain_position)
  {
    d->questionable_associations.erase(
        std::remove_if(d->questionable_associations.begin(),
                       d->questionable_associations.end(),
                       [track_index, chain_position](const InModuleThreadData::QuestionableAssociation& q)
                       {
                         return q.track_index == track_index &&
                                q.chain_position >= chain_position;
                       }),
        d->questionable_associations.end());
  }

  void rollback_track_from_questionable(InModuleThreadData* d,
                                        unsigned int track_index,
                                        unsigned int challenged_blob)
  {
    if (track_index >= d->tracks.size())
    {
      return;
    }

    auto& track = d->tracks[track_index];
    if (!track.active || track.blob_indices.empty())
    {
      return;
    }

    int position = track.questionable_start_position;
    if (position < 0 || static_cast<unsigned int>(position) >= track.blob_indices.size())
    {
      const auto it = std::find(track.blob_indices.begin(), track.blob_indices.end(), challenged_blob);
      if (it == track.blob_indices.end())
      {
        return;
      }
      position = static_cast<int>(std::distance(track.blob_indices.begin(), it));
    }

    for (unsigned int i = static_cast<unsigned int>(position); i < track.blob_indices.size(); ++i)
    {
      auto& bl = d->blobs[track.blob_indices[i]];
      if (bl.owner_tracklet == static_cast<int32_t>(track_index))
      {
        bl.ownership = InModuleThreadData::BLOB_FREE;
        bl.owner_tracklet = -1;
        bl.owner_cost = std::numeric_limits<float>::infinity();
      }
    }

    track.blob_indices.resize(static_cast<unsigned int>(position));
    track.has_questionable = 0;
    track.questionable_start_position = -1;
    track.parked = 1;
    track.needs_repair = 1;
    track.stop_reason = InModuleThreadData::STOP_AMBIGUOUS_NOT_ACCEPTED;
    remove_questionable_records_from(d, track_index, static_cast<unsigned int>(position));

    if (track.blob_indices.size() >= 2)
    {
      InModuleThreadData::Track refit = track;
      if (make_track_from_blob_chain(d->blobs, track.blob_indices,
                                     d->weight_power,
                                     d->adc_weight_floor_frac,
                                     refit))
      {
        refit.track_id = track.track_id;
        refit.pass = track.pass;
        refit.active = 1;
        refit.parked = 1;
        refit.has_questionable = 0;
        refit.needs_repair = 1;
        refit.looper_candidate = track.looper_candidate;
        refit.stop_reason = track.stop_reason;
        refit.questionable_start_position = -1;
        track = std::move(refit);
      }
    }
    else
    {
      track.nblobs = static_cast<unsigned int>(track.blob_indices.size());
      track.nrawhits = count_raw_hits_from_blob_chain(d->blobs, track.blob_indices);
      if (track.blob_indices.empty())
      {
        track.active = 0;
      }
    }
  }

  bool assign_candidate(InModuleThreadData* d,
                        unsigned int track_index,
                        int direction,
                        unsigned int pass,
                        const LayerSearchResult& search,
                        const CandidateEval& selected)
  {
    if (track_index >= d->tracks.size())
    {
      return false;
    }

    auto& track = d->tracks[track_index];
    auto& bl = d->blobs[selected.blob_index];

    if (selected.ownership_relation == 1 && bl.owner_tracklet >= 0)
    {
      const unsigned int old_owner = static_cast<unsigned int>(bl.owner_tracklet);
      rollback_track_from_questionable(d, old_owner, selected.blob_index);
    }

    const bool ambiguous = search.broad_count > 1U;
    const bool already_questionable = track.has_questionable != 0U;
    const bool make_soft = d->enable_questionable_reassignment &&
                           (ambiguous || already_questionable || selected.ownership_relation == 1);

    unsigned int insertion_position = 0;
    if (direction > 0)
    {
      insertion_position = static_cast<unsigned int>(track.blob_indices.size());
      track.blob_indices.push_back(selected.blob_index);
    }
    else
    {
      track.blob_indices.insert(track.blob_indices.begin(), selected.blob_index);
      insertion_position = 0;
      if (track.questionable_start_position >= 0)
      {
        ++track.questionable_start_position;
      }
      for (auto& qa : d->questionable_associations)
      {
        if (qa.track_index == track_index)
        {
          ++qa.chain_position;
        }
      }
    }

    if (selected.ownership_relation == 2)
    {
      bl.ownership = InModuleThreadData::BLOB_SHARED;
      bl.owner_tracklet = -1;
      bl.owner_cost = static_cast<float>(selected.cost);
      track.has_questionable = 1;
      if (track.questionable_start_position < 0)
      {
        track.questionable_start_position = static_cast<int>(insertion_position);
      }
    }
    else if (make_soft)
    {
      bl.ownership = InModuleThreadData::BLOB_SOFT;
      bl.owner_tracklet = static_cast<int32_t>(track_index);
      bl.owner_cost = static_cast<float>(selected.cost);
      if (!track.has_questionable)
      {
        track.questionable_start_position = static_cast<int>(insertion_position);
      }
      track.has_questionable = 1;
    }
    else
    {
      bl.ownership = InModuleThreadData::BLOB_HARD;
      bl.owner_tracklet = static_cast<int32_t>(track_index);
      bl.owner_cost = static_cast<float>(selected.cost);
    }

    if (ambiguous && d->enable_questionable_reassignment)
    {
      InModuleThreadData::QuestionableAssociation qa;
      qa.track_index = track_index;
      qa.blob_index = selected.blob_index;
      qa.chain_position = insertion_position;
      qa.cost = static_cast<float>(selected.cost);
      qa.competitor_cost = static_cast<float>(search.second_cost);
      qa.competitor_blob = -1;
      double best_other = std::numeric_limits<double>::max();
      for (const auto& candidate : search.candidates)
      {
        if (candidate.blob_index == selected.blob_index || !candidate.accepted)
        {
          continue;
        }
        if (candidate.cost < best_other)
        {
          best_other = candidate.cost;
          qa.competitor_blob = static_cast<int>(candidate.blob_index);
        }
      }
      d->questionable_associations.push_back(qa);
    }

    track.pass = static_cast<uint8_t>(std::max<unsigned int>(track.pass, pass));
    track.parked = 0;
    track.stop_reason = InModuleThreadData::STOP_NONE;
    return true;
  }

  void set_track_stop(InModuleThreadData* d,
                      InModuleThreadData::Track& track,
                      uint8_t reason)
  {
    track.parked = 1;
    track.stop_reason = reason;
    if (d->do_pattern_qa)
    {
      d->pattern_stop_reasons.push_back(reason);
    }
  }

  bool grow_one_direction_new(InModuleThreadData* d,
                              unsigned int track_index,
                              int direction,
                              unsigned int pass)
  {
    if (track_index >= d->tracks.size())
    {
      return false;
    }

    bool added_any = false;
    while (true)
    {
      auto& track = d->tracks[track_index];
      if (!track.active || track.blob_indices.empty())
      {
        break;
      }

      const unsigned int edge_layer = direction > 0
                                          ? d->blobs[track.blob_indices.back()].layer
                                          : d->blobs[track.blob_indices.front()].layer;
      const unsigned int first_layer = module_first_layer(d);
      const unsigned int last_layer = module_last_layer(d);
      if ((direction > 0 && edge_layer >= last_layer) ||
          (direction < 0 && edge_layer <= first_layer))
      {
        set_track_stop(d, track, InModuleThreadData::STOP_MODULE_EDGE);
        break;
      }

      bool accepted_step = false;
      bool saw_broad = false;
      for (unsigned int layer_step = 1; layer_step <= d->max_layer_step; ++layer_step)
      {
        const int target_signed = static_cast<int>(edge_layer) +
                                  direction * static_cast<int>(layer_step);
        if (target_signed < static_cast<int>(first_layer) ||
            target_signed > static_cast<int>(last_layer))
        {
          break;
        }
        const unsigned int target_layer = static_cast<unsigned int>(target_signed);

        LayerSearchResult search = search_layer(d, track_index,
                                                d->tracks[track_index].blob_indices,
                                                target_layer, direction, pass);
        record_layer_search_qa(d, search, pass, layer_step);

        if (search.broad_count == 0)
        {
          continue;
        }

        saw_broad = true;
        if (search.best_index < 0)
        {
          auto& current = d->tracks[track_index];
          set_track_stop(d, current,
                         search.broad_count == 1
                             ? InModuleThreadData::STOP_UNIQUE_BROAD_NOT_ACCEPTED
                             : InModuleThreadData::STOP_AMBIGUOUS_NOT_ACCEPTED);
          break;
        }

        const CandidateEval selected = search.candidates[static_cast<unsigned int>(search.best_index)];
        if (!assign_candidate(d, track_index, direction, pass, search, selected))
        {
          auto& current = d->tracks[track_index];
          set_track_stop(d, current, InModuleThreadData::STOP_AMBIGUOUS_NOT_ACCEPTED);
          break;
        }

        accepted_step = true;
        added_any = true;
        break;
      }

      if (!accepted_step)
      {
        auto& current = d->tracks[track_index];
        if (!saw_broad && current.stop_reason == InModuleThreadData::STOP_NONE)
        {
          set_track_stop(d, current, InModuleThreadData::STOP_NO_BROAD);
        }
        break;
      }
    }

    return added_any;
  }

  // -------------------------------------------------------------------
  // Exact legacy path used for the GPFS consistency check. This preserves the
  // original ADC seed order, adjacent-layer growth, one 6x6-like window, and
  // original candidate score. New storage/indexing is transparent to physics.
  // -------------------------------------------------------------------
  int find_best_blob_legacy(const InModuleThreadData* d,
                            unsigned int target_layer,
                            double pred_pad,
                            double pred_tbin)
  {
    int best = -1;
    double best_score = std::numeric_limits<double>::max();
    const auto range = blob_range_for_layer(d, target_layer);

    for (unsigned int i = range.first; i < range.second; ++i)
    {
      const auto& bl = d->blobs[i];
      if (bl.ownership != InModuleThreadData::BLOB_FREE)
      {
        continue;
      }

      const double dp = bl.pad - pred_pad;
      const double dt = bl.tbin - pred_tbin;
      if (std::fabs(dp) > d->pass1_tight_dp ||
          std::fabs(dt) > d->pass1_tight_dt)
      {
        continue;
      }

      const double score = sqr(dp) / (sqr(static_cast<double>(d->pass1_tight_dp)) + 1.0e-9) +
                           sqr(dt) / (sqr(static_cast<double>(d->pass1_tight_dt)) + 1.0e-9) -
                           0.01 * std::log(bl.adc + 1.0);
      if (score < best_score)
      {
        best_score = score;
        best = static_cast<int>(i);
      }
    }
    return best;
  }

  void grow_one_direction_legacy(InModuleThreadData* d,
                                 std::vector<unsigned int>& chain,
                                 int direction)
  {
    while (true)
    {
      unsigned int edge_layer = direction > 0
                                    ? d->blobs[chain.back()].layer
                                    : d->blobs[chain.front()].layer;

      if ((direction > 0 && edge_layer >= module_last_layer(d)) ||
          (direction < 0 && edge_layer <= module_first_layer(d)))
      {
        break;
      }

      const unsigned int target_layer = static_cast<unsigned int>(
          static_cast<int>(edge_layer) + direction);

      double pred_pad = direction > 0
                            ? d->blobs[chain.back()].pad
                            : d->blobs[chain.front()].pad;
      double pred_tbin = direction > 0
                             ? d->blobs[chain.back()].tbin
                             : d->blobs[chain.front()].tbin;

      if (chain.size() >= 2)
      {
        InModuleThreadData::Track tmp;
        if (fit_track_from_blobs(d->blobs, chain,
                                 d->weight_power,
                                 d->adc_weight_floor_frac,
                                 tmp))
        {
          pred_pad = tmp.pad_slope * static_cast<double>(target_layer) + tmp.pad_intercept;
          pred_tbin = tmp.tbin_slope * static_cast<double>(target_layer) + tmp.tbin_intercept;
        }
      }

      const int ibest = find_best_blob_legacy(d, target_layer, pred_pad, pred_tbin);
      if (ibest < 0)
      {
        break;
      }

      auto& bl = d->blobs[static_cast<unsigned int>(ibest)];
      bl.ownership = InModuleThreadData::BLOB_HARD;
      bl.owner_tracklet = -2;
      bl.owner_cost = 0.0F;
      if (direction > 0)
      {
        chain.push_back(static_cast<unsigned int>(ibest));
      }
      else
      {
        chain.insert(chain.begin(), static_cast<unsigned int>(ibest));
      }
    }
  }

  void build_tracks_legacy(InModuleThreadData* d)
  {
    d->tracks.clear();
    for (auto& bl : d->blobs)
    {
      bl.ownership = InModuleThreadData::BLOB_FREE;
      bl.owner_tracklet = -1;
      bl.owner_cost = std::numeric_limits<float>::infinity();
    }

    std::vector<unsigned int> order(d->blobs.size());
    std::iota(order.begin(), order.end(), 0U);
    std::sort(order.begin(), order.end(), BlobAdcSort(&d->blobs));

    unsigned int tid = 0;
    for (unsigned int seed : order)
    {
      if (d->blobs[seed].ownership != InModuleThreadData::BLOB_FREE)
      {
        continue;
      }

      std::vector<unsigned int> chain{seed};
      d->blobs[seed].ownership = InModuleThreadData::BLOB_HARD;
      d->blobs[seed].owner_tracklet = -2;
      d->blobs[seed].owner_cost = 0.0F;

      grow_one_direction_legacy(d, chain, +1);
      grow_one_direction_legacy(d, chain, -1);

      if (chain.size() < d->min_track_blobs)
      {
        for (unsigned int k : chain)
        {
          d->blobs[k].ownership = InModuleThreadData::BLOB_FREE;
          d->blobs[k].owner_tracklet = -1;
          d->blobs[k].owner_cost = std::numeric_limits<float>::infinity();
        }
        continue;
      }

      InModuleThreadData::Track trk;
      trk.track_id = tid;
      if (make_track_from_blob_chain(d->blobs, chain,
                                     d->weight_power,
                                     d->adc_weight_floor_frac,
                                     trk))
      {
        trk.track_id = tid++;
        trk.pass = 1;
        trk.active = 1;
        d->tracks.push_back(std::move(trk));
      }
    }
  }

  bool legacy_mode(const InModuleThreadData* d)
  {
    return !d->enable_second_pass &&
           !d->enable_third_pass &&
           !d->enable_questionable_reassignment &&
           !d->seed_from_inner_layers &&
           d->max_layer_step == 1 &&
           d->pass1_tight_dt == d->pass1_broad_dt &&
           d->pass1_tight_dp == d->pass1_broad_dp &&
           d->pass1_max_dpad_slope_residual < 0.0 &&
           d->pass1_max_dtbin_slope_residual < 0.0 &&
           d->pass1_max_dpad_second_residual < 0.0 &&
           d->pass1_max_dtbin_second_residual < 0.0;
  }

  void seed_track(InModuleThreadData* d,
                  unsigned int blob_index,
                  unsigned int pass)
  {
    InModuleThreadData::Track trk;
    trk.track_id = static_cast<unsigned int>(d->tracks.size());
    trk.pass = static_cast<uint8_t>(pass);
    trk.active = 1;
    trk.parked = 0;
    trk.blob_indices.push_back(blob_index);
    trk.nblobs = 1;
    trk.nrawhits = d->blobs[blob_index].raw_hit_count;

    d->tracks.push_back(std::move(trk));
    const unsigned int track_index = static_cast<unsigned int>(d->tracks.size() - 1);
    auto& bl = d->blobs[blob_index];
    bl.ownership = InModuleThreadData::BLOB_HARD;
    bl.owner_tracklet = static_cast<int32_t>(track_index);
    bl.owner_cost = 0.0F;
  }

  void release_track_ownership(InModuleThreadData* d, unsigned int track_index)
  {
    if (track_index >= d->tracks.size())
    {
      return;
    }
    auto& track = d->tracks[track_index];
    for (unsigned int ib : track.blob_indices)
    {
      auto& bl = d->blobs[ib];
      if (bl.owner_tracklet == static_cast<int32_t>(track_index))
      {
        bl.ownership = InModuleThreadData::BLOB_FREE;
        bl.owner_tracklet = -1;
        bl.owner_cost = std::numeric_limits<float>::infinity();
      }
    }
    track.active = 0;
    track.blob_indices.clear();
  }

  double track_curvature_metric(const InModuleThreadData* d,
                                const InModuleThreadData::Track& track)
  {
    if (track.blob_indices.size() < 3)
    {
      return 0.0;
    }
    const DerivativeState ds = compute_derivative_state(d, track.blob_indices, +1,
                                                        d->blobs[track.blob_indices.back()].layer);
    if (!ds.second_valid)
    {
      return 0.0;
    }
    return std::fabs(ds.mean_pad_second) + std::fabs(ds.mean_tbin_second);
  }

  void grow_existing_tracks_pass2(InModuleThreadData* d)
  {
    std::vector<unsigned int> order;
    for (unsigned int i = 0; i < d->tracks.size(); ++i)
    {
      const auto& track = d->tracks[i];
      if (track.active && track.parked &&
          track.stop_reason != InModuleThreadData::STOP_MODULE_EDGE &&
          !track.blob_indices.empty())
      {
        order.push_back(i);
      }
    }

    std::sort(order.begin(), order.end(),
              [d](unsigned int a, unsigned int b)
              {
                const auto& ta = d->tracks[a];
                const auto& tb = d->tracks[b];
                if (ta.blob_indices.size() != tb.blob_indices.size())
                {
                  return ta.blob_indices.size() > tb.blob_indices.size();
                }
                const double ca = track_curvature_metric(d, ta);
                const double cb = track_curvature_metric(d, tb);
                if (ca != cb)
                {
                  return ca < cb;
                }
                return a < b;
              });

    for (unsigned int i : order)
    {
      if (!d->tracks[i].active)
      {
        continue;
      }
      d->tracks[i].pass = 2;
      d->tracks[i].parked = 0;
      d->tracks[i].stop_reason = InModuleThreadData::STOP_NONE;
      grow_one_direction_new(d, i, +1, 2);
      grow_one_direction_new(d, i, -1, 2);
    }
  }

  void seed_remaining_tracks_pass2(InModuleThreadData* d)
  {
    std::vector<unsigned int> order(d->blobs.size());
    std::iota(order.begin(), order.end(), 0U);
    if (d->seed_from_inner_layers)
    {
      std::sort(order.begin(), order.end(), BlobInnerSort(&d->blobs));
    }
    else
    {
      std::sort(order.begin(), order.end(), BlobAdcSort(&d->blobs));
    }

    for (unsigned int seed : order)
    {
      if (d->blobs[seed].ownership != InModuleThreadData::BLOB_FREE)
      {
        continue;
      }

      seed_track(d, seed, 2);
      const unsigned int track_index = static_cast<unsigned int>(d->tracks.size() - 1);
      grow_one_direction_new(d, track_index, +1, 2);
      grow_one_direction_new(d, track_index, -1, 2);

      if (d->tracks[track_index].blob_indices.size() < 2)
      {
        release_track_ownership(d, track_index);
      }
    }
  }

  void repair_rolled_back_tracks(InModuleThreadData* d)
  {
    for (unsigned int iter = 0; iter < d->questionable_repair_iterations; ++iter)
    {
      bool attempted_any = false;
      for (unsigned int i = 0; i < d->tracks.size(); ++i)
      {
        auto& track = d->tracks[i];
        if (!track.active || !track.needs_repair || track.blob_indices.empty())
        {
          continue;
        }

        attempted_any = true;
        track.parked = 0;
        track.needs_repair = 0;
        track.stop_reason = InModuleThreadData::STOP_NONE;
        grow_one_direction_new(d, i, +1, 2);
        grow_one_direction_new(d, i, -1, 2);
      }

      if (!attempted_any)
      {
        break;
      }

      // A repair can steal another SOFT blob and roll back a track that was
      // already visited in this iteration. Run another iteration only when
      // such a pending local repair remains.
      const bool pending = std::any_of(
          d->tracks.begin(), d->tracks.end(),
          [](const InModuleThreadData::Track& track)
          {
            return track.active && track.needs_repair && !track.blob_indices.empty();
          });
      if (!pending)
      {
        break;
      }
    }
  }

  bool third_pass_extend_back(InModuleThreadData* d, unsigned int track_index)
  {
    if (track_index >= d->tracks.size())
    {
      return false;
    }
    auto& track = d->tracks[track_index];
    if (!track.active || track.blob_indices.size() < d->third_min_blobs)
    {
      return false;
    }

    bool added = false;
    for (unsigned int istep = 0; istep < d->third_max_steps; ++istep)
    {
      if (track.blob_indices.size() < 2)
      {
        break;
      }

      const auto& previous = d->blobs[track.blob_indices[track.blob_indices.size() - 2]];
      const auto& last = d->blobs[track.blob_indices.back()];
      const double pred_pad = last.pad + (last.pad - previous.pad);
      const double pred_tbin = last.tbin + (last.tbin - previous.tbin);

      int best_blob = -1;
      double best_score = std::numeric_limits<double>::max();
      unsigned int broad_count = 0;

      const int min_layer = std::max(static_cast<int>(module_first_layer(d)),
                                     static_cast<int>(last.layer) - static_cast<int>(d->third_layer_span));
      const int max_layer = std::min(static_cast<int>(module_last_layer(d)),
                                     static_cast<int>(last.layer) + static_cast<int>(d->third_layer_span));

      for (int layer = min_layer; layer <= max_layer; ++layer)
      {
        const auto range = blob_range_for_layer(d, static_cast<unsigned int>(layer));
        for (unsigned int ib = range.first; ib < range.second; ++ib)
        {
          if (contains_blob(track.blob_indices, ib))
          {
            continue;
          }
          const auto& bl = d->blobs[ib];
          if (bl.ownership != InModuleThreadData::BLOB_FREE &&
              bl.owner_tracklet != static_cast<int32_t>(track_index))
          {
            continue;
          }

          const double dp = bl.pad - pred_pad;
          const double dt = bl.tbin - pred_tbin;
          if (std::fabs(dp) > d->third_broad_dp ||
              std::fabs(dt) > d->third_broad_dt)
          {
            continue;
          }
          ++broad_count;
          if (std::fabs(dp) > d->third_tight_dp ||
              std::fabs(dt) > d->third_tight_dt)
          {
            continue;
          }

          const double score = sqr(dp / safe_scale(d->third_tight_dp)) +
                               sqr(dt / safe_scale(d->third_tight_dt));
          if (score < best_score)
          {
            best_score = score;
            best_blob = static_cast<int>(ib);
          }
        }
      }

      if (best_blob < 0)
      {
        track.stop_reason = broad_count > 0
                                ? InModuleThreadData::STOP_AMBIGUOUS_NOT_ACCEPTED
                                : InModuleThreadData::STOP_NO_BROAD;
        break;
      }

      auto& bl = d->blobs[static_cast<unsigned int>(best_blob)];
      if (d->enable_questionable_reassignment && broad_count > 1)
      {
        bl.ownership = InModuleThreadData::BLOB_SOFT;
        bl.owner_tracklet = static_cast<int32_t>(track_index);
        bl.owner_cost = static_cast<float>(best_score);
        if (!track.has_questionable)
        {
          track.questionable_start_position = static_cast<int>(track.blob_indices.size());
        }
        track.has_questionable = 1;
      }
      else
      {
        bl.ownership = InModuleThreadData::BLOB_HARD;
        bl.owner_tracklet = static_cast<int32_t>(track_index);
        bl.owner_cost = static_cast<float>(best_score);
      }
      track.blob_indices.push_back(static_cast<unsigned int>(best_blob));
      track.pass = 3;
      track.looper_candidate = 1;
      track.parked = 0;
      added = true;
    }

    if (track.pass == 3 && track.stop_reason == InModuleThreadData::STOP_NONE)
    {
      track.stop_reason = InModuleThreadData::STOP_THIRD_PASS_LIMIT;
    }
    return added;
  }

  void run_third_pass(InModuleThreadData* d)
  {
    for (unsigned int i = 0; i < d->tracks.size(); ++i)
    {
      if (!d->tracks[i].active)
      {
        continue;
      }
      const auto& track = d->tracks[i];
      if (track.stop_reason == InModuleThreadData::STOP_MODULE_EDGE ||
          track.blob_indices.size() < d->third_min_blobs ||
          track_curvature_metric(d, track) < d->third_min_curvature_metric)
      {
        continue;
      }
      third_pass_extend_back(d, i);
    }
  }

  void build_tracks_multi_pass(InModuleThreadData* d)
  {
    d->tracks.clear();
    d->questionable_associations.clear();
    d->pattern_qa.clear();
    d->pattern_stop_reasons.clear();
    d->pattern_curvature_metrics.clear();

    for (auto& bl : d->blobs)
    {
      bl.ownership = InModuleThreadData::BLOB_FREE;
      bl.owner_tracklet = -1;
      bl.owner_cost = std::numeric_limits<float>::infinity();
    }

    std::vector<unsigned int> order(d->blobs.size());
    std::iota(order.begin(), order.end(), 0U);
    if (d->seed_from_inner_layers)
    {
      std::sort(order.begin(), order.end(), BlobInnerSort(&d->blobs));
    }
    else
    {
      std::sort(order.begin(), order.end(), BlobAdcSort(&d->blobs));
    }

    for (unsigned int seed : order)
    {
      if (d->blobs[seed].ownership != InModuleThreadData::BLOB_FREE)
      {
        continue;
      }

      seed_track(d, seed, 1);
      const unsigned int track_index = static_cast<unsigned int>(d->tracks.size() - 1);
      grow_one_direction_new(d, track_index, +1, 1);
      grow_one_direction_new(d, track_index, -1, 1);

      if (d->tracks[track_index].blob_indices.size() < 2)
      {
        release_track_ownership(d, track_index);
      }
      else if (!d->enable_second_pass && !d->enable_third_pass &&
               d->tracks[track_index].blob_indices.size() < d->min_track_blobs)
      {
        release_track_ownership(d, track_index);
      }
    }

    if (d->enable_second_pass)
    {
      grow_existing_tracks_pass2(d);
      seed_remaining_tracks_pass2(d);
      if (d->enable_questionable_reassignment)
      {
        repair_rolled_back_tracks(d);
      }
    }

    if (d->do_pattern_qa)
    {
      for (const auto& track : d->tracks)
      {
        if (track.active && track.blob_indices.size() >= 3 &&
            track.stop_reason != InModuleThreadData::STOP_MODULE_EDGE)
        {
          d->pattern_curvature_metrics.push_back(
              static_cast<float>(track_curvature_metric(d, track)));
        }
      }
    }

    if (d->enable_third_pass)
    {
      run_third_pass(d);
    }

    std::vector<InModuleThreadData::Track> tracklets;
    tracklets.reserve(d->tracks.size());
    for (auto& work : d->tracks)
    {
      if (!work.active || work.blob_indices.size() < d->min_tracklet_blobs_for_connection)
      {
        continue;
      }

      InModuleThreadData::Track trk = work;
      if (!make_track_from_blob_chain(d->blobs, work.blob_indices,
                                      d->weight_power,
                                      d->adc_weight_floor_frac,
                                      trk))
      {
        continue;
      }
      trk.track_id = static_cast<unsigned int>(tracklets.size());
      trk.pass = work.pass;
      trk.has_questionable = work.has_questionable;
      trk.looper_candidate = work.looper_candidate;
      tracklets.push_back(std::move(trk));
    }

    d->tracks.swap(tracklets);
    connect_track_pieces_in_module(d);

    d->tracks.erase(
        std::remove_if(d->tracks.begin(), d->tracks.end(),
                       [d](const InModuleThreadData::Track& trk)
                       {
                         return trk.nblobs < d->min_track_blobs;
                       }),
        d->tracks.end());

    for (unsigned int i = 0; i < d->tracks.size(); ++i)
    {
      d->tracks[i].track_id = i;
    }
  }

  void* ProcessModule(void* arg)
  {
    InModuleThreadData* d = static_cast<InModuleThreadData*>(arg);
    if (!d)
    {
      return nullptr;
    }

    collect_raw_hits(d);
    reject_long_pad_noise(d);
    build_blobs(d);

    if (legacy_mode(d))
    {
      build_tracks_legacy(d);
      connect_track_pieces_in_module(d);
    }
    else
    {
      build_tracks_multi_pass(d);
    }

    if (d->verbosity > 1)
    {
      std::cout << "Tpc_ModuleTrackReco worker: region=" << d->region
                << " sector=" << d->sector << " side=" << d->side
                << " raw_hits=" << d->raw_hits.size()
                << " blobs=" << d->blobs.size()
                << " tracks=" << d->tracks.size() << std::endl;
    }

    return nullptr;
  }

}  // anonymous namespace

// ===================================================================
// Struct constructors
// ===================================================================
InModuleThreadData::LayerHitSet::LayerHitSet()
  : layer(0)
  , hitsetkey(0)
  , hitset(nullptr)
{
}

InModuleThreadData::RawHit::RawHit()
  : layer(0)
  , hitsetkey(0)
  , hitkey(0)
  , pad(0)
  , tbin(0)
  , adc(0)
{
}

InModuleThreadData::Blob::Blob()
  : layer(0)
  , pad(0.0)
  , tbin(0.0)
  , adc(0.0)
  , nhits(0)
  , raw_hit_begin(0)
  , raw_hit_count(0)
  , ownership(BLOB_FREE)
  , owner_tracklet(-1)
  , owner_cost(std::numeric_limits<float>::infinity())
{
}

InModuleThreadData::Track::Track()
  : track_id(0)
  , first_layer(0)
  , last_layer(0)
  , nblobs(0)
  , nrawhits(0)
  , pad_slope(0.0)
  , pad_intercept(0.0)
  , tbin_slope(0.0)
  , tbin_intercept(0.0)
  , pass(0)
  , active(1)
  , parked(0)
  , has_questionable(0)
  , looper_candidate(0)
  , needs_repair(0)
  , stop_reason(STOP_NONE)
  , questionable_start_position(-1)
{
}

InModuleThreadData::InModuleThreadData()
  : region(0)
  , sector(0)
  , side(0)
  , module_key(0)
  , pedestal(74.4)
  , verbosity(0)
  , noise_max_consecutive_timebins(10)
  , noise_keep_first_timebins(3)
  , noise_adc_tolerance(5)
  , blob_dt(2)
  , blob_dp(2)
  , enable_second_pass(false)
  , enable_third_pass(false)
  , enable_questionable_reassignment(false)
  , seed_from_inner_layers(false)
  , do_pattern_qa(false)
  , max_layer_step(1)
  , pass1_tight_dt(6)
  , pass1_tight_dp(6)
  , pass1_broad_dt(6)
  , pass1_broad_dp(6)
  , pass2_broad_dt(12)
  , pass2_broad_dp(12)
  , pass2_straight_dt(6)
  , pass2_straight_dp(6)
  , pass2_local_dt(6)
  , pass2_local_dp(6)
  , pass2_combined_dt(8)
  , pass2_combined_dp(8)
  , pass1_max_dpad_slope_residual(-1.0)
  , pass1_max_dtbin_slope_residual(-1.0)
  , pass1_max_dpad_second_residual(-1.0)
  , pass1_max_dtbin_second_residual(-1.0)
  , pass2_max_dpad_slope_residual(-1.0)
  , pass2_max_dtbin_slope_residual(-1.0)
  , pass2_max_dpad_second_residual(-1.0)
  , pass2_max_dtbin_second_residual(-1.0)
  , derivative_history(3)
  , score_dp_scale(6.0)
  , score_dt_scale(6.0)
  , score_dpad_slope_scale(2.0)
  , score_dtbin_slope_scale(2.0)
  , score_dpad_second_scale(2.0)
  , score_dtbin_second_scale(2.0)
  , transfer_cost_ratio(0.7)
  , share_cost_ratio(1.25)
  , questionable_repair_iterations(2)
  , third_tight_dt(6)
  , third_tight_dp(6)
  , third_broad_dt(10)
  , third_broad_dp(10)
  , third_layer_span(1)
  , third_max_steps(16)
  , third_min_blobs(3)
  , third_min_curvature_metric(0.0)
  , min_track_blobs(4)
  , min_tracklet_blobs_for_connection(4)
  , connect_max_layer_gap(8)
  , connect_dp(8.0)
  , connect_dt(8.0)
  , connect_dpad_slope(2.0)
  , connect_dtbin_slope(2.0)
  , weight_power(0.5)
  , adc_weight_floor_frac(0.15)
{
}

// ===================================================================
// Tpc_ModuleTrackReco
// ===================================================================
Tpc_ModuleTrackReco::Tpc_ModuleTrackReco(const std::string& name,
                                         const std::string& filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_patternQAFileName("Tpc_ModuleTrackRecoPatternQA.root")
  , m_outputFile(nullptr)
  , m_patternQAFile(nullptr)
  , m_tree(nullptr)
  , m_hits(nullptr)
  , m_tpcModuleTrackContainer(nullptr)
  , m_event(0)
  , m_maxThreads(72)
  , m_pedestal(0.0)
  , m_noiseMaxConsecutiveTimebins(10)
  , m_noiseKeepFirstTimebins(3)
  , m_noiseAdcTolerance(5)
  , m_blob_dt(2)
  , m_blob_dp(2)
  , m_enableSecondPass(false)
  , m_enableThirdPass(false)
  , m_enableQuestionableReassignment(false)
  , m_seedFromInnerLayers(false)
  , m_doPatternQA(false)
  , m_maxLayerStep(1)
  , m_pass1Tight_dt(6)
  , m_pass1Tight_dp(6)
  , m_pass1Broad_dt(6)
  , m_pass1Broad_dp(6)
  , m_pass2Broad_dt(12)
  , m_pass2Broad_dp(12)
  , m_pass2Straight_dt(6)
  , m_pass2Straight_dp(6)
  , m_pass2Local_dt(6)
  , m_pass2Local_dp(6)
  , m_pass2Combined_dt(8)
  , m_pass2Combined_dp(8)
  , m_pass1MaxDpadSlopeResidual(-1.0)
  , m_pass1MaxDtbinSlopeResidual(-1.0)
  , m_pass1MaxDpadSecondResidual(-1.0)
  , m_pass1MaxDtbinSecondResidual(-1.0)
  , m_pass2MaxDpadSlopeResidual(-1.0)
  , m_pass2MaxDtbinSlopeResidual(-1.0)
  , m_pass2MaxDpadSecondResidual(-1.0)
  , m_pass2MaxDtbinSecondResidual(-1.0)
  , m_derivativeHistory(3)
  , m_score_dp_scale(6.0)
  , m_score_dt_scale(6.0)
  , m_score_dpad_slope_scale(2.0)
  , m_score_dtbin_slope_scale(2.0)
  , m_score_dpad_second_scale(2.0)
  , m_score_dtbin_second_scale(2.0)
  , m_transferCostRatio(0.7)
  , m_shareCostRatio(1.25)
  , m_questionableRepairIterations(2)
  , m_thirdTight_dt(6)
  , m_thirdTight_dp(6)
  , m_thirdBroad_dt(10)
  , m_thirdBroad_dp(10)
  , m_thirdLayerSpan(1)
  , m_thirdMaxSteps(16)
  , m_thirdMinBlobs(3)
  , m_thirdMinCurvatureMetric(0.0)
  , m_minTrackBlobs(4)
  , m_minTrackletBlobsForConnection(4)
  , m_connectMaxLayerGap(8)
  , m_connect_dp(8.0)
  , m_connect_dt(8.0)
  , m_connect_dpad_slope(2.0)
  , m_connect_dtbin_slope(2.0)
  , m_h_pr_straight_dp(nullptr)
  , m_h_pr_straight_dt(nullptr)
  , m_h_pr_local_dp(nullptr)
  , m_h_pr_local_dt(nullptr)
  , m_h_pr_first_pad_residual(nullptr)
  , m_h_pr_first_tbin_residual(nullptr)
  , m_h_pr_second_pad_residual(nullptr)
  , m_h_pr_second_tbin_residual(nullptr)
  , m_h_pr_candidate_count(nullptr)
  , m_h_pr_layer_step(nullptr)
  , m_h_pr_cost(nullptr)
  , m_h_pr_cost_separation(nullptr)
  , m_h_pr_association_state(nullptr)
  , m_h_pr_stop_reason(nullptr)
  , m_h_pr_curvature_metric(nullptr)
  , m_h_pr_straight_vs_local_pad(nullptr)
  , m_h_pr_straight_vs_local_tbin(nullptr)
  , m_tree_event(0)
{
}

Tpc_ModuleTrackReco::~Tpc_ModuleTrackReco()
{
  if (m_patternQAFile)
  {
    m_patternQAFile->Close();
    delete m_patternQAFile;
    m_patternQAFile = nullptr;
  }
  if (m_outputFile)
  {
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = nullptr;
  }
}

void Tpc_ModuleTrackReco::setMaxThreads(unsigned int n)
{
  m_maxThreads = (n == 0) ? 1 : n;
}

void Tpc_ModuleTrackReco::create_pattern_qa_histograms()
{
  if (!m_patternQAFile)
  {
    return;
  }
  m_patternQAFile->cd();
  m_h_pr_straight_dp = new TH1D("h_pr_straight_dp", "straight-fit pad residual;#Delta pad;candidate tests", 240, -24.0, 24.0);
  m_h_pr_straight_dt = new TH1D("h_pr_straight_dt", "straight-fit tbin residual;#Delta tbin;candidate tests", 240, -24.0, 24.0);
  m_h_pr_local_dp = new TH1D("h_pr_local_dp", "local-curvature pad residual;#Delta pad;candidate tests", 240, -24.0, 24.0);
  m_h_pr_local_dt = new TH1D("h_pr_local_dt", "local-curvature tbin residual;#Delta tbin;candidate tests", 240, -24.0, 24.0);
  m_h_pr_first_pad_residual = new TH1D("h_pr_first_pad_residual", "first-derivative pad residual;#Delta(dpad/dlayer);candidate tests", 240, -20.0, 20.0);
  m_h_pr_first_tbin_residual = new TH1D("h_pr_first_tbin_residual", "first-derivative tbin residual;#Delta(dtbin/dlayer);candidate tests", 240, -20.0, 20.0);
  m_h_pr_second_pad_residual = new TH1D("h_pr_second_pad_residual", "second-derivative pad residual;#Delta(d^{2}pad/dlayer^{2});candidate tests", 240, -20.0, 20.0);
  m_h_pr_second_tbin_residual = new TH1D("h_pr_second_tbin_residual", "second-derivative tbin residual;#Delta(d^{2}tbin/dlayer^{2});candidate tests", 240, -20.0, 20.0);
  m_h_pr_candidate_count = new TH1D("h_pr_candidate_count", "broad candidates on tested layer;N candidates;layer tests", 16, -0.5, 15.5);
  m_h_pr_layer_step = new TH1D("h_pr_layer_step", "tested layer step;|#Delta layer|;candidate tests", 8, 0.5, 8.5);
  m_h_pr_cost = new TH1D("h_pr_cost", "association cost;cost;candidate tests", 240, 0.0, 60.0);
  m_h_pr_cost_separation = new TH1D("h_pr_cost_separation", "runner-up minus best cost;second-best - best;selected candidates", 240, 0.0, 60.0);
  m_h_pr_association_state = new TH1D("h_pr_association_state", "selected association state;state;selected candidates", 5, -0.5, 4.5);
  m_h_pr_stop_reason = new TH1D("h_pr_stop_reason", "track growth stop reason;reason;stops", 8, -0.5, 7.5);
  m_h_pr_curvature_metric = new TH1D("h_pr_curvature_metric", "stopped-track curvature metric;|<d^{2}pad/dlayer^{2}>|+|<d^{2}tbin/dlayer^{2}>|;tracks", 240, 0.0, 24.0);
  m_h_pr_straight_vs_local_pad = new TH2D("h_pr_straight_vs_local_pad", "straight vs local pad residual;straight #Delta pad;local #Delta pad", 160, -16.0, 16.0, 160, -16.0, 16.0);
  m_h_pr_straight_vs_local_tbin = new TH2D("h_pr_straight_vs_local_tbin", "straight vs local tbin residual;straight #Delta tbin;local #Delta tbin", 160, -16.0, 16.0, 160, -16.0, 16.0);
}

void Tpc_ModuleTrackReco::fill_pattern_qa(const InModuleThreadData& td)
{
  if (!m_doPatternQA)
  {
    return;
  }

  for (const auto& q : td.pattern_qa)
  {
    if (m_h_pr_straight_dp) m_h_pr_straight_dp->Fill(q.straight_dp);
    if (m_h_pr_straight_dt) m_h_pr_straight_dt->Fill(q.straight_dt);
    if (q.local_valid && m_h_pr_local_dp) m_h_pr_local_dp->Fill(q.local_dp);
    if (q.local_valid && m_h_pr_local_dt) m_h_pr_local_dt->Fill(q.local_dt);
    if (q.first_valid && m_h_pr_first_pad_residual) m_h_pr_first_pad_residual->Fill(q.first_pad_residual);
    if (q.first_valid && m_h_pr_first_tbin_residual) m_h_pr_first_tbin_residual->Fill(q.first_tbin_residual);
    if (q.second_valid && m_h_pr_second_pad_residual) m_h_pr_second_pad_residual->Fill(q.second_pad_residual);
    if (q.second_valid && m_h_pr_second_tbin_residual) m_h_pr_second_tbin_residual->Fill(q.second_tbin_residual);
    if (m_h_pr_candidate_count) m_h_pr_candidate_count->Fill(q.candidate_count);
    if (m_h_pr_layer_step) m_h_pr_layer_step->Fill(q.layer_step);
    if (m_h_pr_cost && std::isfinite(q.cost)) m_h_pr_cost->Fill(q.cost);
    if (q.selected && m_h_pr_cost_separation && std::isfinite(q.second_cost) && std::isfinite(q.cost))
    {
      m_h_pr_cost_separation->Fill(std::max(0.0F, q.second_cost - q.cost));
    }
    if (q.selected && m_h_pr_association_state) m_h_pr_association_state->Fill(q.association_state);
    if (q.local_valid && m_h_pr_straight_vs_local_pad) m_h_pr_straight_vs_local_pad->Fill(q.straight_dp, q.local_dp);
    if (q.local_valid && m_h_pr_straight_vs_local_tbin) m_h_pr_straight_vs_local_tbin->Fill(q.straight_dt, q.local_dt);
  }

  for (uint8_t reason : td.pattern_stop_reasons)
  {
    if (m_h_pr_stop_reason) m_h_pr_stop_reason->Fill(reason);
  }
  for (float curvature : td.pattern_curvature_metrics)
  {
    if (m_h_pr_curvature_metric) m_h_pr_curvature_metric->Fill(curvature);
  }
}

void Tpc_ModuleTrackReco::write_pattern_qa_histograms()
{
  if (!m_patternQAFile)
  {
    return;
  }
  m_patternQAFile->cd();
  if (m_h_pr_straight_dp) m_h_pr_straight_dp->Write();
  if (m_h_pr_straight_dt) m_h_pr_straight_dt->Write();
  if (m_h_pr_local_dp) m_h_pr_local_dp->Write();
  if (m_h_pr_local_dt) m_h_pr_local_dt->Write();
  if (m_h_pr_first_pad_residual) m_h_pr_first_pad_residual->Write();
  if (m_h_pr_first_tbin_residual) m_h_pr_first_tbin_residual->Write();
  if (m_h_pr_second_pad_residual) m_h_pr_second_pad_residual->Write();
  if (m_h_pr_second_tbin_residual) m_h_pr_second_tbin_residual->Write();
  if (m_h_pr_candidate_count) m_h_pr_candidate_count->Write();
  if (m_h_pr_layer_step) m_h_pr_layer_step->Write();
  if (m_h_pr_cost) m_h_pr_cost->Write();
  if (m_h_pr_cost_separation) m_h_pr_cost_separation->Write();
  if (m_h_pr_association_state) m_h_pr_association_state->Write();
  if (m_h_pr_stop_reason) m_h_pr_stop_reason->Write();
  if (m_h_pr_curvature_metric) m_h_pr_curvature_metric->Write();
  if (m_h_pr_straight_vs_local_pad) m_h_pr_straight_vs_local_pad->Write();
  if (m_h_pr_straight_vs_local_tbin) m_h_pr_straight_vs_local_tbin->Write();
}

int Tpc_ModuleTrackReco::Init(PHCompositeNode* /*unused*/)
{
  if (Verbosity() > 0)
  {
    m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
    if (!m_outputFile || m_outputFile->IsZombie())
    {
      std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_tree = new TTree("Tpc_ModuleTrackReco", "TPC in-module pattern recognition");
    m_tree->Branch("event", &m_tree_event, "event/I");
    m_tree->Branch("track_id", &m_tree_track_id);
    m_tree->Branch("region", &m_tree_region);
    m_tree->Branch("sector", &m_tree_sector);
    m_tree->Branch("side", &m_tree_side);
    m_tree->Branch("nblobs", &m_tree_nblobs);
    m_tree->Branch("nrawhits", &m_tree_nrawhits);
    m_tree->Branch("first_layer", &m_tree_first_layer);
    m_tree->Branch("last_layer", &m_tree_last_layer);
    m_tree->Branch("hit_event", &m_tree_hit_event);
    m_tree->Branch("hit_track_id", &m_tree_hit_track_id);
    m_tree->Branch("hit_region", &m_tree_hit_region);
    m_tree->Branch("hit_sector", &m_tree_hit_sector);
    m_tree->Branch("hit_side", &m_tree_hit_side);
    m_tree->Branch("hit_layer", &m_tree_hit_layer);
    m_tree->Branch("hit_hitsetkey", &m_tree_hit_hitsetkey);
    m_tree->Branch("hit_hitkey", &m_tree_hit_hitkey);

    std::cout << Name() << "::Init - output file " << m_outputFileName << " created" << std::endl;
  }

  if (m_doPatternQA)
  {
    m_patternQAFile = new TFile(m_patternQAFileName.c_str(), "RECREATE");
    if (!m_patternQAFile || m_patternQAFile->IsZombie())
    {
      std::cerr << Name() << "::Init - cannot create pattern QA file "
                << m_patternQAFileName << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    create_pattern_qa_histograms();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_ModuleTrackReco::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_ModuleTrackReco::End(PHCompositeNode* /*unused*/)
{
  if (m_outputFile)
  {
    m_outputFile->cd();
    if (m_tree)
    {
      m_tree->Write();
    }
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = nullptr;
  }

  if (m_patternQAFile)
  {
    write_pattern_qa_histograms();
    m_patternQAFile->Close();
    delete m_patternQAFile;
    m_patternQAFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_ModuleTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_ModuleTrackReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode =
      dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_tpcModuleTrackContainer =
      findNode::getClass<Tpc_ModuleTrackContainer>(topNode, "TPC_MODULETRACKS");

  if (!m_tpcModuleTrackContainer)
  {
    m_tpcModuleTrackContainer = new Tpc_ModuleTrackContainerv1();

    PHIODataNode<PHObject>* node =
        new PHIODataNode<PHObject>(m_tpcModuleTrackContainer,
                                   "TPC_MODULETRACKS", "PHObject");
    dstNode->addNode(node);

    std::cout << Name() << "::createNodes - created TPC_MODULETRACKS node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void Tpc_ModuleTrackReco::reset_tree_vars()
{
  m_tree_event = m_event;

  m_tree_track_id.clear();
  m_tree_region.clear();
  m_tree_sector.clear();
  m_tree_side.clear();
  m_tree_nblobs.clear();
  m_tree_nrawhits.clear();
  m_tree_first_layer.clear();
  m_tree_last_layer.clear();

  m_tree_hit_event.clear();
  m_tree_hit_track_id.clear();
  m_tree_hit_region.clear();
  m_tree_hit_sector.clear();
  m_tree_hit_side.clear();
  m_tree_hit_layer.clear();
  m_tree_hit_hitsetkey.clear();
  m_tree_hit_hitkey.clear();
}

int Tpc_ModuleTrackReco::process_event(PHCompositeNode* /*unused*/)
{
  reset_tree_vars();

  if (m_tpcModuleTrackContainer)
  {
    m_tpcModuleTrackContainer->Reset();
  }

  std::vector<InModuleThreadData> tdata;
  tdata.reserve(72);

  for (unsigned int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int region = 0; region < 3; ++region)
      {
        InModuleThreadData td;
        td.region = region;
        td.sector = sector;
        td.side = static_cast<int>(side);
        td.module_key = TpcDefs::genModuleHitSetKey(static_cast<uint8_t>(region),
                                                    static_cast<uint8_t>(sector),
                                                    static_cast<uint8_t>(side));
        td.pedestal = m_pedestal;
        td.verbosity = Verbosity();
        td.noise_max_consecutive_timebins = m_noiseMaxConsecutiveTimebins;
        td.noise_keep_first_timebins = m_noiseKeepFirstTimebins;
        td.noise_adc_tolerance = m_noiseAdcTolerance;
        td.blob_dt = m_blob_dt;
        td.blob_dp = m_blob_dp;

        td.enable_second_pass = m_enableSecondPass;
        td.enable_third_pass = m_enableThirdPass;
        td.enable_questionable_reassignment = m_enableQuestionableReassignment;
        td.seed_from_inner_layers = m_seedFromInnerLayers;
        td.do_pattern_qa = m_doPatternQA;
        td.max_layer_step = m_maxLayerStep;

        td.pass1_tight_dt = m_pass1Tight_dt;
        td.pass1_tight_dp = m_pass1Tight_dp;
        td.pass1_broad_dt = m_pass1Broad_dt;
        td.pass1_broad_dp = m_pass1Broad_dp;

        td.pass2_broad_dt = m_pass2Broad_dt;
        td.pass2_broad_dp = m_pass2Broad_dp;
        td.pass2_straight_dt = m_pass2Straight_dt;
        td.pass2_straight_dp = m_pass2Straight_dp;
        td.pass2_local_dt = m_pass2Local_dt;
        td.pass2_local_dp = m_pass2Local_dp;
        td.pass2_combined_dt = m_pass2Combined_dt;
        td.pass2_combined_dp = m_pass2Combined_dp;

        td.pass1_max_dpad_slope_residual = m_pass1MaxDpadSlopeResidual;
        td.pass1_max_dtbin_slope_residual = m_pass1MaxDtbinSlopeResidual;
        td.pass1_max_dpad_second_residual = m_pass1MaxDpadSecondResidual;
        td.pass1_max_dtbin_second_residual = m_pass1MaxDtbinSecondResidual;
        td.pass2_max_dpad_slope_residual = m_pass2MaxDpadSlopeResidual;
        td.pass2_max_dtbin_slope_residual = m_pass2MaxDtbinSlopeResidual;
        td.pass2_max_dpad_second_residual = m_pass2MaxDpadSecondResidual;
        td.pass2_max_dtbin_second_residual = m_pass2MaxDtbinSecondResidual;
        td.derivative_history = m_derivativeHistory;

        td.score_dp_scale = m_score_dp_scale;
        td.score_dt_scale = m_score_dt_scale;
        td.score_dpad_slope_scale = m_score_dpad_slope_scale;
        td.score_dtbin_slope_scale = m_score_dtbin_slope_scale;
        td.score_dpad_second_scale = m_score_dpad_second_scale;
        td.score_dtbin_second_scale = m_score_dtbin_second_scale;
        td.transfer_cost_ratio = m_transferCostRatio;
        td.share_cost_ratio = m_shareCostRatio;
        td.questionable_repair_iterations = m_questionableRepairIterations;

        td.third_tight_dt = m_thirdTight_dt;
        td.third_tight_dp = m_thirdTight_dp;
        td.third_broad_dt = m_thirdBroad_dt;
        td.third_broad_dp = m_thirdBroad_dp;
        td.third_layer_span = m_thirdLayerSpan;
        td.third_max_steps = m_thirdMaxSteps;
        td.third_min_blobs = m_thirdMinBlobs;
        td.third_min_curvature_metric = m_thirdMinCurvatureMetric;

        td.min_track_blobs = m_minTrackBlobs;
        td.min_tracklet_blobs_for_connection = m_minTrackletBlobsForConnection;
        td.connect_max_layer_gap = m_connectMaxLayerGap;
        td.connect_dp = m_connect_dp;
        td.connect_dt = m_connect_dt;
        td.connect_dpad_slope = m_connect_dpad_slope;
        td.connect_dtbin_slope = m_connect_dtbin_slope;

        for (unsigned int l = 0; l < 16; ++l)
        {
          const unsigned int layer = region * 16 + l + 7;
          const TrkrDefs::hitsetkey hitset_key = TpcDefs::genHitSetKey(layer, sector, side);
          TrkrHitSet* hitset = m_hits->findHitSet(hitset_key);
          if (!hitset)
          {
            continue;
          }

          InModuleThreadData::LayerHitSet lhs;
          lhs.layer = layer;
          lhs.hitsetkey = hitset_key;
          lhs.hitset = hitset;
          td.layer_hitsets.push_back(lhs);
        }

        if (!td.layer_hitsets.empty())
        {
          tdata.push_back(std::move(td));
        }
      }
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " has " << tdata.size() << " non-empty modules" << std::endl;
  }

  const unsigned int maxLive = std::max(1U,
                                        std::min(m_maxThreads,
                                                 static_cast<unsigned int>(tdata.size())));

  for (unsigned int start = 0; start < static_cast<unsigned int>(tdata.size()); start += maxLive)
  {
    const unsigned int end = std::min(start + maxLive,
                                      static_cast<unsigned int>(tdata.size()));
    const unsigned int nLive = end - start;

    std::vector<pthread_t> threads(nLive);
    std::vector<int> thread_ok(nLive, 0);

    for (unsigned int i = 0; i < nLive; ++i)
    {
      const unsigned int idx = start + i;
      const int rc = pthread_create(&threads[i], nullptr, ProcessModule,
                                    static_cast<void*>(&tdata[idx]));
      if (rc != 0)
      {
        std::cerr << Name() << "::process_event - pthread_create failed for"
                  << " region=" << tdata[idx].region
                  << " sector=" << tdata[idx].sector
                  << " side=" << tdata[idx].side << std::endl;
      }
      else
      {
        thread_ok[i] = 1;
      }
    }

    for (unsigned int i = 0; i < nLive; ++i)
    {
      if (thread_ok[i])
      {
        pthread_join(threads[i], nullptr);
      }
    }
  }

  // Harvest results and QA only after every worker has joined.
  for (const auto& td : tdata)
  {
    fill_pattern_qa(td);

    for (const auto& tr : td.tracks)
    {
      const unsigned int global_track_id =
          m_tpcModuleTrackContainer
              ? m_tpcModuleTrackContainer->size()
              : static_cast<unsigned int>(m_tree_track_id.size());

      Tpc_ModuleTrackv1* outTrack = new Tpc_ModuleTrackv1();
      outTrack->set_event(static_cast<unsigned int>(m_event));
      outTrack->set_track_id(global_track_id);
      outTrack->set_region(td.region);
      outTrack->set_sector(td.sector);
      outTrack->set_side(td.side);
      outTrack->set_nblobs(tr.nblobs);
      outTrack->set_nrawhits(tr.nrawhits);
      outTrack->set_first_layer(tr.first_layer);
      outTrack->set_last_layer(tr.last_layer);

      m_tree_track_id.push_back(global_track_id);
      m_tree_region.push_back(td.region);
      m_tree_sector.push_back(td.sector);
      m_tree_side.push_back(td.side);
      m_tree_nblobs.push_back(tr.nblobs);
      m_tree_nrawhits.push_back(tr.nrawhits);
      m_tree_first_layer.push_back(tr.first_layer);
      m_tree_last_layer.push_back(tr.last_layer);

      for (unsigned int ib : tr.blob_indices)
      {
        const auto& bl = td.blobs[ib];
        const unsigned int end = bl.raw_hit_begin + bl.raw_hit_count;
        for (unsigned int iflat = bl.raw_hit_begin; iflat < end; ++iflat)
        {
          const unsigned int raw_hit_index = td.blob_raw_hit_indices[iflat];
          const InModuleThreadData::RawHit& rh = td.raw_hits[raw_hit_index];

          outTrack->add_hit_index(rh.hitsetkey, rh.hitkey);

          m_tree_hit_event.push_back(static_cast<unsigned int>(m_event));
          m_tree_hit_track_id.push_back(global_track_id);
          m_tree_hit_region.push_back(td.region);
          m_tree_hit_sector.push_back(td.sector);
          m_tree_hit_side.push_back(td.side);
          m_tree_hit_layer.push_back(rh.layer);
          m_tree_hit_hitsetkey.push_back(static_cast<unsigned long long>(rh.hitsetkey));
          m_tree_hit_hitkey.push_back(static_cast<unsigned long long>(rh.hitkey));
        }
      }

      if (m_tpcModuleTrackContainer)
      {
        m_tpcModuleTrackContainer->add_track(outTrack);
      }
      else
      {
        delete outTrack;
      }
    }
  }

  if (m_tree)
  {
    m_tree->Fill();
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " tracks=" << m_tree_track_id.size()
              << " track-raw-hits=" << m_tree_hit_track_id.size() << std::endl;
  }
  if (m_tpcModuleTrackContainer && Verbosity() > 0)
  {
    m_tpcModuleTrackContainer->identify();
  }
  if (m_tpcModuleTrackContainer && Verbosity() > 1)
  {
    for (unsigned int i = 0; i < m_tpcModuleTrackContainer->size(); ++i)
    {
      const Tpc_ModuleTrack* trk = m_tpcModuleTrackContainer->get_track(i);
      if (trk)
      {
        trk->identify();
      }
    }
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
