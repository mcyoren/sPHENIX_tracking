// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCASSEMBLEDTRACKRECO_H
#define TPCTRACKRECO_TPCASSEMBLEDTRACKRECO_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

class IdealPadMap;

#include <cstdint>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class TH1D;
class TH2D;

class Tpc_ModuleTrackContainer;
class Tpc_AssembledTrackContainer;
class TrkrHitSetContainer;

class Tpc_AssembledTrackReco : public SubsysReco
{
 public:
  explicit Tpc_AssembledTrackReco(const std::string& name = "Tpc_AssembledTrackReco",
                                  const std::string& filename = "Tpc_AssembledTracks.root");
  ~Tpc_AssembledTrackReco() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setDebugOutputFileName(const std::string& n) { m_debugOutputFileName = n; }
  void setDoDebugHistograms(bool v) { m_doDebugHistograms = v; }

  // Applied only after full assembly, before writing TPC_ASSEMBLEDTRACKS.
  void setMinClustersPerTrack(unsigned int n) { m_minClustersPerTrack = n; }

  // Loose linear-phi preselection before exact sagitta matching.
  void setLinearPhiPrecutScale(double v) { m_linearPhiPrecutScale = v; }

  void setConnectMaxLayerGap(unsigned int n) { m_connectMaxLayerGap = n; }

  // Connection window in global phi radians and tbin units at match radius.
  void setConnectWindow(double dphi, double dtbin)
  {
    m_connect_dphi = dphi;
    m_connect_dtbin = dtbin;
  }

  // Slope windows: d(phi)/d(radius) and d(tbin)/d(radius).
  void setConnectSlopeWindow(double dphi_slope, double dtbin_slope)
  {
    m_connect_dphi_slope = dphi_slope;
    m_connect_dtbin_slope = dtbin_slope;
  }

  void setUseSagittaPhiFit(bool v) { m_useSagittaPhiFit = v; }
  void setUseAnalyticSagittaSlope(bool v) { m_useAnalyticSagittaSlope = v; }

  // Optional fast stage. Defaults are OFF so the old assembly path remains
  // available for exact consistency tests.
  void setFastSameSectorPass(bool v) { m_enableFastSameSectorPass = v; }
  void setFastSameSectorScales(double broad_scale, double tight_scale)
  {
    m_fastBroadScale = broad_scale;
    m_fastTightScale = tight_scale;
  }

  // Optional unique-tight shortcut in the global normal stage. When enabled,
  // one unique cheap candidate satisfying the tight linear/tbin windows can be
  // trial-refit without first doing the exact pairwise sagitta test.
  void setNormalUniqueTightPass(bool v) { m_enableNormalUniqueTight = v; }
  void setNormalTightScale(double v) { m_normalTightScale = v; }

  // Optional rescue. Long-long merges remain protected; a short candidate can
  // still be attached to a long candidate. Scale grows linearly over iterations.
  void setRescuePass(bool enabled,
                     unsigned int iterations = 1,
                     double max_window_scale = 1.5,
                     unsigned int max_segments = 2)
  {
    m_enableRescuePass = enabled;
    m_rescueIterations = iterations;
    m_rescueMaxWindowScale = max_window_scale;
    m_rescueMaxSegments = max_segments;
  }

  void setSeedCovarianceDiagonal(double sigmaX, double sigmaY, double sigmaZ,
                                 double sigmaPx, double sigmaPy, double sigmaPz)
  {
    m_seedSigmaX = sigmaX;
    m_seedSigmaY = sigmaY;
    m_seedSigmaZ = sigmaZ;
    m_seedSigmaPx = sigmaPx;
    m_seedSigmaPy = sigmaPy;
    m_seedSigmaPz = sigmaPz;
  }

 public:
  struct PiecePoint
  {
    double radius{0.0};
    double phi{0.0};
    double tbin{0.0};
    double weight{0.0};

    TrkrDefs::hitsetkey hitsetkey{0};
    TrkrDefs::hitkey hitkey{0};
  };

  struct Piece
  {
    Piece();

    unsigned int source_index;
    unsigned int source_track_id;
    unsigned int event;
    unsigned int region;
    unsigned int sector;
    int side;

    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int nblobs;
    unsigned int nrawhits;
    uint64_t layer_mask;

    double phi_slope;
    double phi_intercept;
    mutable double phi_S;
    mutable double phi_x0;
    mutable double phi_invR;
    mutable double phi_theta;
    mutable double phi_bline;
    mutable bool phi_sagitta_ok;
    mutable bool phi_sagitta_evaluated;

    double tbin_slope;
    double tbin_intercept;

    std::vector<PiecePoint> points;
  };

  struct SeedParameters
  {
    bool ok{false};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double px{0.0};
    double py{0.0};
    double pz{0.0};
    double cov[6][6]{};
  };

  struct Candidate
  {
    Candidate();

    unsigned int event;
    int side;
    unsigned int first_layer;
    unsigned int last_layer;
    unsigned int first_sector;
    unsigned int last_sector;
    unsigned int first_region;
    unsigned int last_region;
    unsigned int nsegments;
    unsigned int nblobs;
    unsigned int nrawhits;
    uint64_t layer_mask;

    double phi_slope;
    double phi_intercept;
    double phi_S;
    double phi_x0;
    double phi_invR;
    double phi_theta;
    double phi_bline;
    bool phi_sagitta_ok;

    double tbin_slope_r;
    double tbin_intercept_r;
    double chi2_phi;
    double chi2_tbin;
    int ndof_phi;
    int ndof_tbin;

    std::vector<unsigned int> piece_indices;
  };

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void reset_tree_vars();
  void create_debug_histograms();
  void write_debug_histograms();

  bool make_piece(unsigned int source_index, Piece& p) const;
  void ensure_piece_sagitta(const Piece& p) const;

  double predict_phi(const Piece& p, double radius) const;
  double predict_phi(const Candidate& c, double radius) const;
  double predict_phi_slope(const Piece& p, double radius) const;
  double predict_phi_slope(const Candidate& c, double radius) const;
  std::pair<double, double> predict_phi_and_slope(const Piece& p, double radius) const;
  std::pair<double, double> predict_phi_and_slope(const Candidate& c, double radius) const;

  bool refit_candidate(const std::vector<Piece>& pieces,
                       const std::vector<unsigned int>& piece_indices,
                       Candidate& c) const;

  bool candidates_can_connect(const Candidate& a,
                              const Piece& b,
                              double& score,
                              double& b_phi_intercept_shifted) const;

  bool candidates_can_connect(const Candidate& a,
                              const Candidate& b,
                              double& score,
                              bool allow_same_sector = false,
                              double window_scale = 1.0) const;

  bool cheap_candidate_relation(const Candidate& a,
                                const Candidate& b,
                                bool allow_same_sector,
                                double window_scale,
                                bool& broad,
                                bool& tight,
                                double& score) const;

  bool fast_piece_relation(const Candidate& a,
                           const Piece& b,
                           bool& broad,
                           bool& tight,
                           double& score) const;

  void connect_sector_pieces(const std::vector<Piece>& pieces,
                             int side,
                             unsigned int sector,
                             std::vector<Candidate>& output) const;

  void connect_sector_pieces_fast(const std::vector<Piece>& pieces,
                                  int side,
                                  unsigned int sector,
                                  std::vector<Candidate>& output) const;

  void connect_side_candidates(const std::vector<Piece>& pieces,
                               const std::vector<Candidate>& seeds,
                               int side,
                               std::vector<Candidate>& output,
                               bool allow_same_sector = false,
                               double window_scale = 1.0,
                               bool enable_unique_tight = false) const;

  void rescue_side_candidates(const std::vector<Piece>& pieces,
                              const std::vector<Candidate>& seeds,
                              int side,
                              double window_scale,
                              std::vector<Candidate>& output) const;

  SeedParameters make_seed_parameters(const Candidate& c) const;

  std::string m_outputFileName;
  std::string m_debugOutputFileName;
  std::string m_inputNodeName;
  std::string m_outputNodeName;

  TFile* m_outputFile;
  TFile* m_debugOutputFile;
  TTree* m_tree;

  Tpc_ModuleTrackContainer* m_tpcModuleTrackContainer;
  Tpc_AssembledTrackContainer* m_assembledTrackContainer;
  TrkrHitSetContainer* m_hits;

  int m_event;
  IdealPadMap* m_idealPadMap;

  unsigned int m_connectMaxLayerGap;
  double m_connect_dphi;
  double m_connect_dtbin;
  double m_connect_dphi_slope;
  double m_connect_dtbin_slope;
  bool m_useSagittaPhiFit;
  bool m_useAnalyticSagittaSlope;
  bool m_doDebugHistograms;
  unsigned int m_minClustersPerTrack;
  double m_linearPhiPrecutScale;

  bool m_enableFastSameSectorPass;
  double m_fastBroadScale;
  double m_fastTightScale;
  bool m_enableNormalUniqueTight;
  double m_normalTightScale;

  bool m_enableRescuePass;
  unsigned int m_rescueIterations;
  double m_rescueMaxWindowScale;
  unsigned int m_rescueMaxSegments;

  double m_seedSigmaX;
  double m_seedSigmaY;
  double m_seedSigmaZ;
  double m_seedSigmaPx;
  double m_seedSigmaPy;
  double m_seedSigmaPz;

  // Debug matching histograms
  TH1D* m_h_dphi;
  TH1D* m_h_dtbin;
  TH1D* m_h_dmphi;
  TH1D* m_h_dmtbin;
  TH1D* m_h_score;

  TH2D* m_h_dphi_vs_dtbin;
  TH2D* m_h_dmphi_vs_dmtbin;
  TH2D* m_h_dphi_vs_dmphi;

  TH2D* m_h_tbin_slope_vs_first_tbin;
  TH2D* m_h_tbin_slope_vs_last_tbin;

  TH2D* m_h_track_tbin_slope_vs_tbin_span_3modules;
  TH2D* m_h_track_tbin_slope_vs_first_tbin_3modules;
  TH2D* m_h_track_tbin_slope_vs_last_tbin_3modules;

  TH1D* m_h_layer_gap;
  TH1D* m_h_nsegments;
  TH1D* m_h_matched_sector_delta;

  mutable std::mutex m_debugMutex;

  int m_tree_event{};

  std::vector<unsigned int> m_tree_track_id;
  std::vector<int> m_tree_side;

  std::vector<unsigned int> m_tree_nsegments;
  std::vector<unsigned int> m_tree_nblobs;
  std::vector<unsigned int> m_tree_nrawhits;

  std::vector<unsigned int> m_tree_first_layer;
  std::vector<unsigned int> m_tree_last_layer;

  std::vector<unsigned int> m_tree_first_sector;
  std::vector<unsigned int> m_tree_last_sector;

  std::vector<unsigned int> m_tree_first_region;
  std::vector<unsigned int> m_tree_last_region;

  std::vector<unsigned int> m_tree_source_assembled_track_id;
  std::vector<unsigned int> m_tree_source_inmodule_track_id;
  std::vector<unsigned int> m_tree_source_region;
  std::vector<unsigned int> m_tree_source_sector;
  std::vector<int> m_tree_source_side;

  std::vector<unsigned int> m_tree_hit_assembled_track_id;
  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};

#endif
