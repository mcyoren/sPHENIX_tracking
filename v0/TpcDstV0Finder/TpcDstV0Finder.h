#pragma once

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class FinalTrack;
class FinalTrackContainer;
class FinalTrackVertexContainer;

class TpcDstV0Finder : public SubsysReco
{
 public:
  explicit TpcDstV0Finder(const std::string& name = "TpcDstV0Finder",
                          const std::string& output = "TpcDstV0.root");
  ~TpcDstV0Finder() override = default;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void set_output_file(const std::string& value) { m_output_file = value; }
  void set_track_node(const std::string& value) { m_track_node = value; }
  void set_vertex_node(const std::string& value) { m_vertex_node = value; }
  void set_bfield(double value) { m_bfield_t = value; }
  void set_primary_vertex(double x, double y, double z);
  void use_vertex_node(bool value = true) { m_use_vertex_node = value; }
  void write_track_tree(bool value = true) { m_write_track_tree = value; }
  void write_same_sign_pairs(bool value = true) { m_write_same_sign = value; }
  void enable_vertex_constrained_refit(bool value = true) { m_refit_after_vertex = value; }

  void set_min_clusters(unsigned int value) { m_min_clusters = value; }
  void set_min_pt(double value) { m_min_pt = value; }
  void set_max_track_chi2_ndf(double value) { m_max_track_chi2_ndf = value; }
  void set_max_pair_dca(double value) { m_max_pair_dca = value; }
  void set_min_decay_radius(double value) { m_min_decay_radius = value; }
  void set_max_abs_decay_z(double value) { m_max_abs_decay_z = value; }
  void set_min_dira(double value) { m_min_dira = value; }
  void set_min_lproj(double value) { m_min_lproj = value; }
  void set_max_abs_alpha(double value) { m_max_abs_alpha = value; }
  void set_search_range(double backward_cm, double forward_cm);
  void set_coarse_steps(unsigned int value) { m_coarse_steps = value; }

 private:
  struct Vec3
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  struct State
  {
    unsigned int id{0};
    unsigned int source_id{0};
    unsigned int nclusters{0};
    int fit_status{0};
    int charge{0};
    double chi2{0.0};
    double ndf{0.0};
    double dedx{0.0};
    Vec3 r0;
    Vec3 p0;
    double pt{0.0};
    double p{0.0};
    double phi0{0.0};
    double tanl{0.0};
    double omega{0.0}; // dphi/ds, s is transverse arc length [cm]
  };

  struct PcaResult
  {
    bool valid{false};
    double s1{0.0};
    double s2{0.0};
    Vec3 x1;
    Vec3 x2;
    double dca{0.0};
  };

  struct PairRow
  {
    int run{0};
    int evt{0};
    unsigned int track_id1{0};
    unsigned int track_id2{0};
    unsigned int source_id1{0};
    unsigned int source_id2{0};
    short charge1{0};
    short charge2{0};
    unsigned short nclusters1{0};
    unsigned short nclusters2{0};
    float quality1{0};
    float quality2{0};
    float dedx_1{0};
    float dedx_2{0};
    float ref_x1{0}, ref_y1{0}, ref_z1{0};
    float ref_x2{0}, ref_y2{0}, ref_z2{0};
    float ref_px1{0}, ref_py1{0}, ref_pz1{0};
    float ref_px2{0}, ref_py2{0}, ref_pz2{0};
    float pca1_x{0}, pca1_y{0}, pca1_z{0};
    float pca2_x{0}, pca2_y{0}, pca2_z{0};
    float pca_x{0}, pca_y{0}, pca_z{0};
    float px1{0}, py1{0}, pz1{0};
    float px2{0}, py2{0}, pz2{0};
    float v0_px{0}, v0_py{0}, v0_pz{0}, v0_pt{0};
    float pairDCA{0};
    float dca_xy1{0}, dca_z1{0}, dca_xy2{0}, dca_z2{0};
    float decay_radius{0};
    float Lproj{0};
    float cosThetaReco{0};
    float alpha{0};
    float qT{0};
    float mass_Kshort{0};
    float mass_Lambda{0};
    float mass_AntiLambda{0};
    float primary_x{0}, primary_y{0}, primary_z{0};
    float s1{0}, s2{0};
    short refit_used{0};
  };

  struct TrackRow
  {
    int run{0};
    int evt{0};
    unsigned int track_id{0};
    unsigned int source_id{0};
    unsigned int nclusters{0};
    int fit_status{0};
    short charge{0};
    float x{0}, y{0}, z{0};
    float px{0}, py{0}, pz{0}, pt{0}, p{0}, eta{0};
    float chi2{0}, ndf{0}, quality{0}, dedx{0};
    float dca_xy{0}, dca_z{0};
    float primary_x{0}, primary_y{0}, primary_z{0};
  };

  bool make_state(const FinalTrack*, State&) const;
  Vec3 position(const State&, double s) const;
  Vec3 momentum(const State&, double s) const;
  Vec3 tangent(const State&, double s) const;
  PcaResult pair_pca(const State&, const State&) const;
  PcaResult refine_pair(const State&, const State&, double s1, double s2) const;
  std::pair<double, double> dca_to_vertex(const State&, const Vec3&) const;
  Vec3 choose_primary_vertex(const FinalTrackVertexContainer*) const;
  bool fill_pair(const State&, const State&, const Vec3&, int run, int evt);
  void fill_track(const State&, const Vec3&, int run, int evt);
  void create_branches();
  int get_run(PHCompositeNode*) const;
  int get_event(PHCompositeNode*) const;

  static Vec3 add(const Vec3&, const Vec3&);
  static Vec3 sub(const Vec3&, const Vec3&);
  static Vec3 scale(const Vec3&, double);
  static double dot(const Vec3&, const Vec3&);
  static double norm(const Vec3&);
  static double distance(const Vec3&, const Vec3&);
  static double cosine(const Vec3&, const Vec3&);
  static double invariant_mass(const Vec3&, double, const Vec3&, double);
  static bool armenteros(const Vec3& pplus, const Vec3& pminus, double& alpha, double& qt);
  static float f(double value);

  std::string m_output_file;
  std::string m_track_node{"FINALTRACKS"};
  std::string m_vertex_node{"FINALTRACKVERTICES"};
  double m_bfield_t{1.4};
  Vec3 m_fixed_primary{0.0, 0.0, 0.0};
  bool m_use_vertex_node{true};
  bool m_write_track_tree{true};
  bool m_write_same_sign{false};
  bool m_refit_after_vertex{false};

  unsigned int m_min_clusters{20};
  double m_min_pt{0.15};
  double m_max_track_chi2_ndf{-1.0};
  double m_max_pair_dca{6.0};
  double m_min_decay_radius{0.0};
  double m_max_abs_decay_z{110.0};
  double m_min_dira{-1.0};
  double m_min_lproj{0.0};
  double m_max_abs_alpha{1.1};
  double m_search_backward_cm{120.0};
  double m_search_forward_cm{120.0};
  unsigned int m_coarse_steps{40};

  TFile* m_file{nullptr};
  TTree* m_pair_tree{nullptr};
  TTree* m_track_tree{nullptr};
  PairRow m_pair;
  TrackRow m_track;

  std::uint64_t m_events{0};
  std::uint64_t m_tracks{0};
  std::uint64_t m_pairs_tested{0};
  std::uint64_t m_pairs_written{0};
};
