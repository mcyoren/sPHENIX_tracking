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

  void setBlobWindow(int dt, int dp)
  {
    m_blob_dt = dt;
    m_blob_dp = dp;
  }

  // Backward-compatible setter. It makes Pass-1 TIGHT and BROAD identical.
  // With second/third pass off, derivative cuts disabled and max layer step=1,
  // this reproduces the old 6x6 adjacent-layer growth logic.
  void setSearchWindow(int dt, int dp)
  {
    m_pass1Tight_dt = dt;
    m_pass1Tight_dp = dp;
    m_pass1Broad_dt = dt;
    m_pass1Broad_dp = dp;
  }

  void setPass1Window(int tight_dt, int tight_dp, int broad_dt, int broad_dp)
  {
    m_pass1Tight_dt = tight_dt;
    m_pass1Tight_dp = tight_dp;
    m_pass1Broad_dt = broad_dt;
    m_pass1Broad_dp = broad_dp;
  }

  void setSecondPass(bool value) { m_enableSecondPass = value; }
  void setThirdPass(bool value) { m_enableThirdPass = value; }
  void setQuestionableReassignment(bool value) { m_enableQuestionableReassignment = value; }
  void setSeedFromInnerLayers(bool value) { m_seedFromInnerLayers = value; }

  // 1 = adjacent layer only, 2 = allow one missing layer, 3 = allow two.
  void setMaxLayerStep(unsigned int value)
  {
    m_maxLayerStep = value > 0 ? value : 1;
  }

  void setMaxMissingLayers(unsigned int value)
  {
    m_maxLayerStep = value + 1;
  }

  void setPass2BroadWindow(int dt, int dp)
  {
    m_pass2Broad_dt = dt;
    m_pass2Broad_dp = dp;
  }

  void setPass2StraightWindow(int dt, int dp)
  {
    m_pass2Straight_dt = dt;
    m_pass2Straight_dp = dp;
  }

  void setPass2LocalWindow(int dt, int dp)
  {
    m_pass2Local_dt = dt;
    m_pass2Local_dp = dp;
  }

  void setPass2CombinedWindow(int dt, int dp)
  {
    m_pass2Combined_dt = dt;
    m_pass2Combined_dp = dp;
  }

  // pass must be 1 or 2. Negative values disable a derivative cut.
  void setDerivativeCuts(unsigned int pass,
                         double dtbin_first,
                         double pad_first,
                         double dtbin_second,
                         double pad_second)
  {
    if (pass == 1)
    {
      m_pass1MaxDtbinSlopeResidual = dtbin_first;
      m_pass1MaxDpadSlopeResidual = pad_first;
      m_pass1MaxDtbinSecondResidual = dtbin_second;
      m_pass1MaxDpadSecondResidual = pad_second;
    }
    else if (pass == 2)
    {
      m_pass2MaxDtbinSlopeResidual = dtbin_first;
      m_pass2MaxDpadSlopeResidual = pad_first;
      m_pass2MaxDtbinSecondResidual = dtbin_second;
      m_pass2MaxDpadSecondResidual = pad_second;
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
    m_score_dt_scale = dt;
    m_score_dp_scale = dp;
    m_score_dtbin_slope_scale = dtbin_first;
    m_score_dpad_slope_scale = pad_first;
    m_score_dtbin_second_scale = dtbin_second;
    m_score_dpad_second_scale = pad_second;
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
    m_thirdTight_dt = tight_dt;
    m_thirdTight_dp = tight_dp;
    m_thirdBroad_dt = broad_dt;
    m_thirdBroad_dp = broad_dp;
  }

  void setThirdPassLayerSpan(unsigned int n) { m_thirdLayerSpan = n; }
  void setThirdPassMaxSteps(unsigned int n) { m_thirdMaxSteps = n; }
  void setThirdPassMinBlobs(unsigned int n) { m_thirdMinBlobs = n; }
  void setThirdPassMinCurvature(double v) { m_thirdMinCurvatureMetric = v; }

  void setDoPatternQA(bool value) { m_doPatternQA = value; }
  void setPatternQAFileName(const std::string& value) { m_patternQAFileName = value; }

  // Reject long same-pad tails before blob/track building.
  void setNoiseRejection(int max_consecutive_timebins = 10,
                         int keep_first_timebins = 3,
                         int adc_tolerance = 5)
  {
    m_noiseMaxConsecutiveTimebins = max_consecutive_timebins;
    m_noiseKeepFirstTimebins = keep_first_timebins;
    m_noiseAdcTolerance = adc_tolerance;
  }

  void setMinTrackBlobs(unsigned int n) { m_minTrackBlobs = n; }
  void setMinTrackletBlobsForConnection(unsigned int n) { m_minTrackletBlobsForConnection = n; }

  void setConnectMaxLayerGap(unsigned int n) { m_connectMaxLayerGap = n; }

  void setConnectWindow(double dt, double dp)
  {
    m_connect_dt = dt;
    m_connect_dp = dp;
  }

  void setConnectSlopeWindow(double dtbin_slope, double dpad_slope)
  {
    m_connect_dtbin_slope = dtbin_slope;
    m_connect_dpad_slope = dpad_slope;
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

  // Noise rejection
  int m_noiseMaxConsecutiveTimebins;
  int m_noiseKeepFirstTimebins;
  int m_noiseAdcTolerance;

  // Blob building
  int m_blob_dt;
  int m_blob_dp;

  // Pattern-recognition configuration. Defaults preserve old behavior.
  bool m_enableSecondPass;
  bool m_enableThirdPass;
  bool m_enableQuestionableReassignment;
  bool m_seedFromInnerLayers;
  bool m_doPatternQA;

  unsigned int m_maxLayerStep;

  int m_pass1Tight_dt;
  int m_pass1Tight_dp;
  int m_pass1Broad_dt;
  int m_pass1Broad_dp;

  int m_pass2Broad_dt;
  int m_pass2Broad_dp;
  int m_pass2Straight_dt;
  int m_pass2Straight_dp;
  int m_pass2Local_dt;
  int m_pass2Local_dp;
  int m_pass2Combined_dt;
  int m_pass2Combined_dp;

  double m_pass1MaxDpadSlopeResidual;
  double m_pass1MaxDtbinSlopeResidual;
  double m_pass1MaxDpadSecondResidual;
  double m_pass1MaxDtbinSecondResidual;

  double m_pass2MaxDpadSlopeResidual;
  double m_pass2MaxDtbinSlopeResidual;
  double m_pass2MaxDpadSecondResidual;
  double m_pass2MaxDtbinSecondResidual;

  unsigned int m_derivativeHistory;

  double m_score_dp_scale;
  double m_score_dt_scale;
  double m_score_dpad_slope_scale;
  double m_score_dtbin_slope_scale;
  double m_score_dpad_second_scale;
  double m_score_dtbin_second_scale;

  double m_transferCostRatio;
  double m_shareCostRatio;
  unsigned int m_questionableRepairIterations;

  int m_thirdTight_dt;
  int m_thirdTight_dp;
  int m_thirdBroad_dt;
  int m_thirdBroad_dp;
  unsigned int m_thirdLayerSpan;
  unsigned int m_thirdMaxSteps;
  unsigned int m_thirdMinBlobs;
  double m_thirdMinCurvatureMetric;

  unsigned int m_minTrackBlobs;
  unsigned int m_minTrackletBlobsForConnection;

  // Track-piece connection parameters
  unsigned int m_connectMaxLayerGap;

  double m_connect_dp;
  double m_connect_dt;

  double m_connect_dpad_slope;
  double m_connect_dtbin_slope;

  // Pattern QA. These are filled only after module worker threads join.
  TH1D* m_h_pr_straight_dp;
  TH1D* m_h_pr_straight_dt;
  TH1D* m_h_pr_local_dp;
  TH1D* m_h_pr_local_dt;
  TH1D* m_h_pr_first_pad_residual;
  TH1D* m_h_pr_first_tbin_residual;
  TH1D* m_h_pr_second_pad_residual;
  TH1D* m_h_pr_second_tbin_residual;
  TH1D* m_h_pr_candidate_count;
  TH1D* m_h_pr_layer_step;
  TH1D* m_h_pr_cost;
  TH1D* m_h_pr_cost_separation;
  TH1D* m_h_pr_association_state;
  TH1D* m_h_pr_stop_reason;
  TH1D* m_h_pr_curvature_metric;
  TH2D* m_h_pr_straight_vs_local_pad;
  TH2D* m_h_pr_straight_vs_local_tbin;

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
