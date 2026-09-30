// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCMODULETRACKRECO_H
#define TPCTRACKRECO_TPCMODULETRACKRECO_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;

class TFile;
class TTree;
class TH1D;
class TH2D;

class Tpc_ModuleTrackContainer;

class TrkrHitSetContainer;
class TrkrHitSet;

// ===================================================================
// Per-module thread data
// ===================================================================
struct InModuleThreadData
{
  InModuleThreadData();

  enum BlobOwnership : uint8_t
  {
    BLOB_FREE = 0,
    BLOB_HARD = 1,
    BLOB_SOFT = 2,
    BLOB_SHARED = 3
  };

  enum StopReason : uint8_t
  {
    STOP_NONE = 0,
    STOP_NO_BROAD = 1,
    STOP_UNIQUE_BROAD_NOT_ACCEPTED = 2,
    STOP_AMBIGUOUS_NOT_ACCEPTED = 3,
    STOP_MODULE_EDGE = 4,
    STOP_THIRD_PASS_LIMIT = 5
  };

  struct LayerHitSet
  {
    LayerHitSet();

    unsigned int layer;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrHitSet* hitset;
  };

  struct RawHit
  {
    RawHit();

    unsigned int layer;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrDefs::hitkey hitkey;

    unsigned short pad;
    unsigned short tbin;
    unsigned short adc;
  };

  struct Blob
  {
    Blob();

    unsigned int layer;

    double pad;
    double tbin;
    double adc;

    unsigned int nhits;

    // Raw-hit membership is stored in one flat per-module array.
    unsigned int raw_hit_begin;
    unsigned int raw_hit_count;

    // Pattern-recognition ownership. SOFT blobs remain available to later
    // tracks; HARD blobs do not.
    uint8_t ownership;
    int32_t owner_tracklet;
    float owner_cost;
  };

  struct Track
  {
    Track();

    unsigned int track_id;

    unsigned int first_layer;
    unsigned int last_layer;

    unsigned int nblobs;
    unsigned int nrawhits;

    // Internal temporary straight-line parameters used only for chain growth
    // and piece connection. They are not copied to Tpc_ModuleTrackContainer or
    // saved as track fit output.
    double pad_slope;
    double pad_intercept;
    double tbin_slope;
    double tbin_intercept;

    uint8_t pass;
    uint8_t active;
    uint8_t parked;
    uint8_t has_questionable;
    uint8_t looper_candidate;
    uint8_t needs_repair;
    uint8_t stop_reason;

    // First questionable position in blob_indices, or -1 if the full chain is
    // hard/confident. If a later track wins a SOFT blob, rollback begins here.
    int questionable_start_position;

    std::vector<unsigned int> blob_indices;
  };

  struct QuestionableAssociation
  {
    unsigned int track_index{0};
    unsigned int blob_index{0};
    unsigned int chain_position{0};
    int competitor_blob{-1};
    float cost{0.0F};
    float competitor_cost{0.0F};
  };

  struct PatternQAEntry
  {
    uint8_t pass{0};
    uint8_t selected{0};
    uint8_t association_state{0};
    uint8_t layer_step{0};
    uint8_t candidate_count{0};
    uint8_t local_valid{0};
    uint8_t first_valid{0};
    uint8_t second_valid{0};

    float straight_dp{0.0F};
    float straight_dt{0.0F};
    float local_dp{0.0F};
    float local_dt{0.0F};

    float first_pad_residual{0.0F};
    float first_tbin_residual{0.0F};
    float second_pad_residual{0.0F};
    float second_tbin_residual{0.0F};

    float cost{0.0F};
    float second_cost{0.0F};
  };

  // Detector/module identity
  unsigned int region;
  unsigned int sector;
  int side;
  TrkrDefs::hitsetkey module_key;

  // General configuration
  double pedestal;
  int verbosity;

  // Noise rejection
  int noise_max_consecutive_timebins;
  int noise_keep_first_timebins;
  int noise_adc_tolerance;

  // Blob building
  int blob_dt;
  int blob_dp;

  // Pattern-recognition switches
  bool enable_second_pass;
  bool enable_third_pass;
  bool enable_questionable_reassignment;
  bool seed_from_inner_layers;
  bool do_pattern_qa;

  // 1 = adjacent layer only (old behavior), 2 = skip one missing broad-empty
  // layer, 3 = skip two missing broad-empty layers.
  unsigned int max_layer_step;

  // Pass 1. BROAD is candidate discovery/ambiguity counting; TIGHT is
  // acceptance.
  int pass1_tight_dt;
  int pass1_tight_dp;
  int pass1_broad_dt;
  int pass1_broad_dp;

  // Pass 2. Candidate discovery is the union around the global straight
  // prediction and the local-curvature prediction.
  int pass2_broad_dt;
  int pass2_broad_dp;
  int pass2_straight_dt;
  int pass2_straight_dp;
  int pass2_local_dt;
  int pass2_local_dp;
  int pass2_combined_dt;
  int pass2_combined_dp;

  // Negative values disable the corresponding derivative requirement.
  double pass1_max_dpad_slope_residual;
  double pass1_max_dtbin_slope_residual;
  double pass1_max_dpad_second_residual;
  double pass1_max_dtbin_second_residual;

  double pass2_max_dpad_slope_residual;
  double pass2_max_dtbin_slope_residual;
  double pass2_max_dpad_second_residual;
  double pass2_max_dtbin_second_residual;

  unsigned int derivative_history;

  // Fixed scales used by the pass-independent association cost. They also
  // make questionable ownership comparisons meaningful across passes.
  double score_dp_scale;
  double score_dt_scale;
  double score_dpad_slope_scale;
  double score_dtbin_slope_scale;
  double score_dpad_second_scale;
  double score_dtbin_second_scale;

  // new_cost < transfer_cost_ratio * old_cost transfers a SOFT association.
  // Similar costs with max/min <= share_cost_ratio can be shared.
  double transfer_cost_ratio;
  double share_cost_ratio;
  unsigned int questionable_repair_iterations;

  // Experimental stage 3: follow already established very-low-pT pieces in
  // sequence space rather than requiring monotonically increasing layer.
  int third_tight_dt;
  int third_tight_dp;
  int third_broad_dt;
  int third_broad_dp;
  unsigned int third_layer_span;
  unsigned int third_max_steps;
  unsigned int third_min_blobs;
  double third_min_curvature_metric;

  unsigned int min_track_blobs;
  unsigned int min_tracklet_blobs_for_connection;

  // Track-piece connection.
  unsigned int connect_max_layer_gap;

  double connect_dp;
  double connect_dt;

  double connect_dpad_slope;
  double connect_dtbin_slope;

  // ADC weighting for temporary internal fits
  double weight_power;
  double adc_weight_floor_frac;

  // Per-module containers
  std::vector<LayerHitSet> layer_hitsets;
  std::vector<RawHit> raw_hits;

  // Flat raw-hit membership for every Blob.
  std::vector<unsigned int> blob_raw_hit_indices;

  // Blobs remain in layer order. For local module layer i=0..15, the range is
  // [blob_layer_offsets[i], blob_layer_offsets[i+1]).
  std::vector<Blob> blobs;
  std::array<unsigned int, 17> blob_layer_offsets{};

  std::vector<Track> tracks;
  std::vector<QuestionableAssociation> questionable_associations;

  // Filled only when do_pattern_qa is true. ROOT is touched only by the main
  // thread after worker joins.
  std::vector<PatternQAEntry> pattern_qa;
  std::vector<uint8_t> pattern_stop_reasons;
  std::vector<float> pattern_curvature_metrics;
};

// ===================================================================
// Main Fun4All module
// ===================================================================
class Tpc_ModuleTrackReco : public SubsysReco
{
 public:
  explicit Tpc_ModuleTrackReco(const std::string& name = "Tpc_ModuleTrackReco",
                               const std::string& filename = "Tpc_ModuleTrackReco.root");

  ~Tpc_ModuleTrackReco() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setMaxThreads(unsigned int n);

  void setPedestal(double p)
  {
    m_pedestal = p;
  }

  // TPC region numbering follows TpcDefs: 0=R1, 1=R2, 2=R3.
  static constexpr unsigned int kRegionCount = 3;

  // Existing setters remain backward compatible: they apply the same value to
  // all three TPC regions. The setRegion* variants allow independent R1/R2/R3
  // tuning without changing reconstruction logic.
  void setBlobWindow(int dt, int dp)
  {
    m_blob_dt.fill(dt);
    m_blob_dp.fill(dp);
  }

  void setRegionBlobWindow(unsigned int region, int dt, int dp)
  {
    if (region >= kRegionCount) return;
    m_blob_dt[region] = dt;
    m_blob_dp[region] = dp;
  }

  // Backward-compatible setter. It makes Pass-1 TIGHT and BROAD identical.
  // With second/third pass off, derivative cuts disabled and max layer step=1,
  // this reproduces the old 6x6 adjacent-layer growth logic.
  void setSearchWindow(int dt, int dp)
  {
    m_pass1Tight_dt.fill(dt);
    m_pass1Tight_dp.fill(dp);
    m_pass1Broad_dt.fill(dt);
    m_pass1Broad_dp.fill(dp);
  }

  void setRegionSearchWindow(unsigned int region, int dt, int dp)
  {
    if (region >= kRegionCount) return;
    m_pass1Tight_dt[region] = dt;
    m_pass1Tight_dp[region] = dp;
    m_pass1Broad_dt[region] = dt;
    m_pass1Broad_dp[region] = dp;
  }

  void setPass1Window(int tight_dt, int tight_dp, int broad_dt, int broad_dp)
  {
    m_pass1Tight_dt.fill(tight_dt);
    m_pass1Tight_dp.fill(tight_dp);
    m_pass1Broad_dt.fill(broad_dt);
    m_pass1Broad_dp.fill(broad_dp);
  }

  void setRegionPass1Window(unsigned int region,
                            int tight_dt, int tight_dp,
                            int broad_dt, int broad_dp)
  {
    if (region >= kRegionCount) return;
    m_pass1Tight_dt[region] = tight_dt;
    m_pass1Tight_dp[region] = tight_dp;
    m_pass1Broad_dt[region] = broad_dt;
    m_pass1Broad_dp[region] = broad_dp;
  }

  void setSecondPass(bool value) { m_enableSecondPass = value; }
  void setThirdPass(bool value) { m_enableThirdPass = value; }
  void setQuestionableReassignment(bool value) { m_enableQuestionableReassignment = value; }
  void setSeedFromInnerLayers(bool value) { m_seedFromInnerLayers = value; }

  // 1 = adjacent layer only, 2 = allow one missing layer, 3 = allow two.
  void setMaxLayerStep(unsigned int value)
  {
    m_maxLayerStep.fill(value > 0 ? value : 1);
  }

  void setRegionMaxLayerStep(unsigned int region, unsigned int value)
  {
    if (region >= kRegionCount) return;
    m_maxLayerStep[region] = value > 0 ? value : 1;
  }

  void setMaxMissingLayers(unsigned int value)
  {
    m_maxLayerStep.fill(value + 1);
  }

  void setRegionMaxMissingLayers(unsigned int region, unsigned int value)
  {
    if (region >= kRegionCount) return;
    m_maxLayerStep[region] = value + 1;
  }

  void setPass2BroadWindow(int dt, int dp)
  {
    m_pass2Broad_dt.fill(dt);
    m_pass2Broad_dp.fill(dp);
  }

  void setRegionPass2BroadWindow(unsigned int region, int dt, int dp)
  {
    if (region >= kRegionCount) return;
    m_pass2Broad_dt[region] = dt;
    m_pass2Broad_dp[region] = dp;
  }

  void setPass2StraightWindow(int dt, int dp)
  {
    m_pass2Straight_dt.fill(dt);
    m_pass2Straight_dp.fill(dp);
  }

  void setRegionPass2StraightWindow(unsigned int region, int dt, int dp)
  {
    if (region >= kRegionCount) return;
    m_pass2Straight_dt[region] = dt;
    m_pass2Straight_dp[region] = dp;
  }

  void setPass2LocalWindow(int dt, int dp)
  {
    m_pass2Local_dt.fill(dt);
    m_pass2Local_dp.fill(dp);
  }

  void setRegionPass2LocalWindow(unsigned int region, int dt, int dp)
  {
    if (region >= kRegionCount) return;
    m_pass2Local_dt[region] = dt;
    m_pass2Local_dp[region] = dp;
  }

  void setPass2CombinedWindow(int dt, int dp)
  {
    m_pass2Combined_dt.fill(dt);
    m_pass2Combined_dp.fill(dp);
  }

  void setRegionPass2CombinedWindow(unsigned int region, int dt, int dp)
  {
    if (region >= kRegionCount) return;
    m_pass2Combined_dt[region] = dt;
    m_pass2Combined_dp[region] = dp;
  }

  // pass must be 1 or 2. Negative values disable a derivative cut.
  void setDerivativeCuts(unsigned int pass,
                         double dtbin_first,
                         double pad_first,
                         double dtbin_second,
                         double pad_second)
  {
    for (unsigned int region = 0; region < kRegionCount; ++region)
    {
      setRegionDerivativeCuts(region, pass, dtbin_first, pad_first,
                              dtbin_second, pad_second);
    }
  }

  void setRegionDerivativeCuts(unsigned int region,
                               unsigned int pass,
                               double dtbin_first,
                               double pad_first,
                               double dtbin_second,
                               double pad_second)
  {
    if (region >= kRegionCount) return;
    if (pass == 1)
    {
      m_pass1MaxDtbinSlopeResidual[region] = dtbin_first;
      m_pass1MaxDpadSlopeResidual[region] = pad_first;
      m_pass1MaxDtbinSecondResidual[region] = dtbin_second;
      m_pass1MaxDpadSecondResidual[region] = pad_second;
    }
    else if (pass == 2)
    {
      m_pass2MaxDtbinSlopeResidual[region] = dtbin_first;
      m_pass2MaxDpadSlopeResidual[region] = pad_first;
      m_pass2MaxDtbinSecondResidual[region] = dtbin_second;
      m_pass2MaxDpadSecondResidual[region] = pad_second;
    }
  }

  void setDerivativeHistory(unsigned int n)
  {
    m_derivativeHistory = n > 0 ? n : 1;
  }

  void setAssociationScoreScales(double dt,
                                 double dp,
                                 double dtbin_first,
                                 double pad_first,
                                 double dtbin_second,
                                 double pad_second)
  {
    for (unsigned int region = 0; region < kRegionCount; ++region)
    {
      setRegionAssociationScoreScales(region, dt, dp, dtbin_first, pad_first,
                                      dtbin_second, pad_second);
    }
  }

  void setRegionAssociationScoreScales(unsigned int region,
                                       double dt,
                                       double dp,
                                       double dtbin_first,
                                       double pad_first,
                                       double dtbin_second,
                                       double pad_second)
  {
    if (region >= kRegionCount) return;
    m_score_dt_scale[region] = dt;
    m_score_dp_scale[region] = dp;
    m_score_dtbin_slope_scale[region] = dtbin_first;
    m_score_dpad_slope_scale[region] = pad_first;
    m_score_dtbin_second_scale[region] = dtbin_second;
    m_score_dpad_second_scale[region] = pad_second;
  }

  void setQuestionableOwnership(double transfer_cost_ratio,
                                double share_cost_ratio)
  {
    m_transferCostRatio = transfer_cost_ratio;
    m_shareCostRatio = share_cost_ratio;
  }

  void setQuestionableRepairIterations(unsigned int n)
  {
    m_questionableRepairIterations = n;
  }

  void setThirdPassWindow(int tight_dt, int tight_dp,
                          int broad_dt, int broad_dp)
  {
    m_thirdTight_dt.fill(tight_dt);
    m_thirdTight_dp.fill(tight_dp);
    m_thirdBroad_dt.fill(broad_dt);
    m_thirdBroad_dp.fill(broad_dp);
  }

  void setRegionThirdPassWindow(unsigned int region,
                                int tight_dt, int tight_dp,
                                int broad_dt, int broad_dp)
  {
    if (region >= kRegionCount) return;
    m_thirdTight_dt[region] = tight_dt;
    m_thirdTight_dp[region] = tight_dp;
    m_thirdBroad_dt[region] = broad_dt;
    m_thirdBroad_dp[region] = broad_dp;
  }

  void setThirdPassLayerSpan(unsigned int n) { m_thirdLayerSpan.fill(n); }
  void setThirdPassMaxSteps(unsigned int n) { m_thirdMaxSteps.fill(n); }
  void setThirdPassMinBlobs(unsigned int n) { m_thirdMinBlobs.fill(n); }
  void setThirdPassMinCurvature(double v) { m_thirdMinCurvatureMetric.fill(v); }

  void setRegionThirdPassLayerSpan(unsigned int region, unsigned int n)
  {
    if (region < kRegionCount) m_thirdLayerSpan[region] = n;
  }
  void setRegionThirdPassMaxSteps(unsigned int region, unsigned int n)
  {
    if (region < kRegionCount) m_thirdMaxSteps[region] = n;
  }
  void setRegionThirdPassMinBlobs(unsigned int region, unsigned int n)
  {
    if (region < kRegionCount) m_thirdMinBlobs[region] = n;
  }
  void setRegionThirdPassMinCurvature(unsigned int region, double v)
  {
    if (region < kRegionCount) m_thirdMinCurvatureMetric[region] = v;
  }

  void setDoPatternQA(bool value) { m_doPatternQA = value; }
  void setPatternQAFileName(const std::string& value) { m_patternQAFileName = value; }

  // Reject long same-pad tails before blob/track building.
  void setNoiseRejection(int max_consecutive_timebins = 10,
                         int keep_first_timebins = 3,
                         int adc_tolerance = 5)
  {
    m_noiseMaxConsecutiveTimebins.fill(max_consecutive_timebins);
    m_noiseKeepFirstTimebins.fill(keep_first_timebins);
    m_noiseAdcTolerance.fill(adc_tolerance);
  }

  void setRegionNoiseRejection(unsigned int region,
                               int max_consecutive_timebins = 10,
                               int keep_first_timebins = 3,
                               int adc_tolerance = 5)
  {
    if (region >= kRegionCount) return;
    m_noiseMaxConsecutiveTimebins[region] = max_consecutive_timebins;
    m_noiseKeepFirstTimebins[region] = keep_first_timebins;
    m_noiseAdcTolerance[region] = adc_tolerance;
  }

  void setMinTrackBlobs(unsigned int n) { m_minTrackBlobs.fill(n); }
  void setMinTrackletBlobsForConnection(unsigned int n) { m_minTrackletBlobsForConnection.fill(n); }

  void setRegionMinTrackBlobs(unsigned int region, unsigned int n)
  {
    if (region < kRegionCount) m_minTrackBlobs[region] = n;
  }
  void setRegionMinTrackletBlobsForConnection(unsigned int region, unsigned int n)
  {
    if (region < kRegionCount) m_minTrackletBlobsForConnection[region] = n;
  }

  void setConnectMaxLayerGap(unsigned int n) { m_connectMaxLayerGap.fill(n); }
  void setRegionConnectMaxLayerGap(unsigned int region, unsigned int n)
  {
    if (region < kRegionCount) m_connectMaxLayerGap[region] = n;
  }

  void setConnectWindow(double dt, double dp)
  {
    m_connect_dt.fill(dt);
    m_connect_dp.fill(dp);
  }

  void setRegionConnectWindow(unsigned int region, double dt, double dp)
  {
    if (region >= kRegionCount) return;
    m_connect_dt[region] = dt;
    m_connect_dp[region] = dp;
  }

  void setConnectSlopeWindow(double dtbin_slope, double dpad_slope)
  {
    m_connect_dtbin_slope.fill(dtbin_slope);
    m_connect_dpad_slope.fill(dpad_slope);
  }

  void setRegionConnectSlopeWindow(unsigned int region,
                                   double dtbin_slope,
                                   double dpad_slope)
  {
    if (region >= kRegionCount) return;
    m_connect_dtbin_slope[region] = dtbin_slope;
    m_connect_dpad_slope[region] = dpad_slope;
  }
 private:
  int getNodes(PHCompositeNode*);
  void reset_tree_vars();
  int createNodes(PHCompositeNode*);

  void create_pattern_qa_histograms();
  void fill_pattern_qa(const InModuleThreadData& td);
  void write_pattern_qa_histograms();

  std::string m_outputFileName;
  std::string m_patternQAFileName;

  TFile* m_outputFile;
  TFile* m_patternQAFile;
  TTree* m_tree;

  TrkrHitSetContainer* m_hits;

  Tpc_ModuleTrackContainer* m_tpcModuleTrackContainer;
  int m_event;
  unsigned int m_maxThreads;

  // General configuration
  double m_pedestal;

  // Region-dependent tuning. Index 0/1/2 corresponds to R1/R2/R3.
  std::array<int, kRegionCount> m_noiseMaxConsecutiveTimebins;
  std::array<int, kRegionCount> m_noiseKeepFirstTimebins;
  std::array<int, kRegionCount> m_noiseAdcTolerance;

  std::array<int, kRegionCount> m_blob_dt;
  std::array<int, kRegionCount> m_blob_dp;

  // Pattern-recognition switches remain common across regions. Numerical cuts
  // and windows below are independent for R1/R2/R3.
  bool m_enableSecondPass;
  bool m_enableThirdPass;
  bool m_enableQuestionableReassignment;
  bool m_seedFromInnerLayers;
  bool m_doPatternQA;

  std::array<unsigned int, kRegionCount> m_maxLayerStep;

  std::array<int, kRegionCount> m_pass1Tight_dt;
  std::array<int, kRegionCount> m_pass1Tight_dp;
  std::array<int, kRegionCount> m_pass1Broad_dt;
  std::array<int, kRegionCount> m_pass1Broad_dp;

  std::array<int, kRegionCount> m_pass2Broad_dt;
  std::array<int, kRegionCount> m_pass2Broad_dp;
  std::array<int, kRegionCount> m_pass2Straight_dt;
  std::array<int, kRegionCount> m_pass2Straight_dp;
  std::array<int, kRegionCount> m_pass2Local_dt;
  std::array<int, kRegionCount> m_pass2Local_dp;
  std::array<int, kRegionCount> m_pass2Combined_dt;
  std::array<int, kRegionCount> m_pass2Combined_dp;

  std::array<double, kRegionCount> m_pass1MaxDpadSlopeResidual;
  std::array<double, kRegionCount> m_pass1MaxDtbinSlopeResidual;
  std::array<double, kRegionCount> m_pass1MaxDpadSecondResidual;
  std::array<double, kRegionCount> m_pass1MaxDtbinSecondResidual;

  std::array<double, kRegionCount> m_pass2MaxDpadSlopeResidual;
  std::array<double, kRegionCount> m_pass2MaxDtbinSlopeResidual;
  std::array<double, kRegionCount> m_pass2MaxDpadSecondResidual;
  std::array<double, kRegionCount> m_pass2MaxDtbinSecondResidual;

  unsigned int m_derivativeHistory;

  std::array<double, kRegionCount> m_score_dp_scale;
  std::array<double, kRegionCount> m_score_dt_scale;
  std::array<double, kRegionCount> m_score_dpad_slope_scale;
  std::array<double, kRegionCount> m_score_dtbin_slope_scale;
  std::array<double, kRegionCount> m_score_dpad_second_scale;
  std::array<double, kRegionCount> m_score_dtbin_second_scale;

  double m_transferCostRatio;
  double m_shareCostRatio;
  unsigned int m_questionableRepairIterations;

  std::array<int, kRegionCount> m_thirdTight_dt;
  std::array<int, kRegionCount> m_thirdTight_dp;
  std::array<int, kRegionCount> m_thirdBroad_dt;
  std::array<int, kRegionCount> m_thirdBroad_dp;
  std::array<unsigned int, kRegionCount> m_thirdLayerSpan;
  std::array<unsigned int, kRegionCount> m_thirdMaxSteps;
  std::array<unsigned int, kRegionCount> m_thirdMinBlobs;
  std::array<double, kRegionCount> m_thirdMinCurvatureMetric;

  std::array<unsigned int, kRegionCount> m_minTrackBlobs;
  std::array<unsigned int, kRegionCount> m_minTrackletBlobsForConnection;

  // Track-piece connection parameters
  std::array<unsigned int, kRegionCount> m_connectMaxLayerGap;
  std::array<double, kRegionCount> m_connect_dp;
  std::array<double, kRegionCount> m_connect_dt;
  std::array<double, kRegionCount> m_connect_dpad_slope;
  std::array<double, kRegionCount> m_connect_dtbin_slope;

  // Pattern QA is stored independently for R1/R2/R3 so the numerical cuts
  // can be tuned from detector-region-specific residual distributions.
  std::array<TH1D*, kRegionCount> m_h_pr_straight_dp{};
  std::array<TH1D*, kRegionCount> m_h_pr_straight_dt{};
  std::array<TH1D*, kRegionCount> m_h_pr_local_dp{};
  std::array<TH1D*, kRegionCount> m_h_pr_local_dt{};
  std::array<TH1D*, kRegionCount> m_h_pr_first_pad_residual{};
  std::array<TH1D*, kRegionCount> m_h_pr_first_tbin_residual{};
  std::array<TH1D*, kRegionCount> m_h_pr_second_pad_residual{};
  std::array<TH1D*, kRegionCount> m_h_pr_second_tbin_residual{};
  std::array<TH1D*, kRegionCount> m_h_pr_candidate_count{};
  std::array<TH1D*, kRegionCount> m_h_pr_layer_step{};
  std::array<TH1D*, kRegionCount> m_h_pr_cost{};
  std::array<TH1D*, kRegionCount> m_h_pr_cost_separation{};
  std::array<TH1D*, kRegionCount> m_h_pr_association_state{};
  std::array<TH1D*, kRegionCount> m_h_pr_stop_reason{};
  std::array<TH1D*, kRegionCount> m_h_pr_curvature_metric{};
  std::array<TH2D*, kRegionCount> m_h_pr_straight_vs_local_pad{};
  std::array<TH2D*, kRegionCount> m_h_pr_straight_vs_local_tbin{};

  // Event number saved once per tree entry
  int m_tree_event;

  // One entry per found module-track. No fit branches are saved here.
  std::vector<unsigned int> m_tree_track_id;
  std::vector<unsigned int> m_tree_region;
  std::vector<unsigned int> m_tree_sector;
  std::vector<int> m_tree_side;

  std::vector<unsigned int> m_tree_nblobs;
  std::vector<unsigned int> m_tree_nrawhits;

  std::vector<unsigned int> m_tree_first_layer;
  std::vector<unsigned int> m_tree_last_layer;

  // Flat per-hit content for TTree reading.
  // Hits are identified by their TrkrHitSetContainer keys only;
  // no hit data is duplicated here.
  std::vector<unsigned int> m_tree_hit_event;
  std::vector<unsigned int> m_tree_hit_track_id;

  std::vector<unsigned int> m_tree_hit_region;
  std::vector<unsigned int> m_tree_hit_sector;
  std::vector<int> m_tree_hit_side;

  std::vector<unsigned int> m_tree_hit_layer;

  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};
#endif
