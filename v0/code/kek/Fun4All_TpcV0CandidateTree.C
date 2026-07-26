#ifndef MACRO_FUN4ALLTPCV0CANDIDATETREE_C
#define MACRO_FUN4ALLTPCV0CANDIDATETREE_C

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)

#include <trackingdiagnostics/TpcV0CandidateTree.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <string>

namespace
{
  std::string env_string(const char *name, const std::string &fallback)
  {
    const char *value = std::getenv(name);
    return value ? std::string(value) : fallback;
  }

  int env_int(const char *name, const int fallback)
  {
    const char *value = std::getenv(name);
    return value ? std::stoi(value) : fallback;
  }

  double env_double(const char *name, const double fallback)
  {
    const char *value = std::getenv(name);
    return value ? std::stod(value) : fallback;
  }

  bool env_bool(const char *name, const bool fallback)
  {
    const std::string value = env_string(name, fallback ? "1" : "0");
    return value == "1" || value == "true" || value == "TRUE" ||
           value == "yes" || value == "YES" || value == "on" || value == "ON";
  }
}

int Fun4All_TpcV0CandidateTree(
    const std::string &input_list = "input.list",
    const std::string &output_file = "TpcV0Candidates.root",
    const int nEvents = 0,
    const int skipEvents = 0)
{

  auto *se = Fun4AllServer::instance();
  se->Verbosity(env_int("V0_SERVER_VERBOSITY", 0));

  auto *v0 = new TpcV0CandidateTree("TpcV0CandidateTree", output_file);
  v0->Verbosity(env_int("V0_MODULE_VERBOSITY", 0));

  v0->use_pattern_cluster_tracks(env_bool("V0_USE_PATTERN_TRACKS", true));
  v0->set_tpc_sa_cluster_node(env_string("V0_CLUSTER_NODE", "TPC_POLYCLUSTERS"));
  v0->set_tpc_sa_track_node(env_string("V0_TRACK_NODE", "TPC_POLYTRACKS"));
  v0->set_tpc_sa_track_vertex_node(env_string("V0_VERTEX_NODE", "TPC_POLYTRACKVERTICES"));

  // Real data: use the fixed beam position, not G4 truth.
  v0->set_use_truth_primary_vertex(env_bool("V0_USE_TRUTH_VERTEX", false));
  v0->set_primary_vertex(
      env_double("V0_BEAM_X", 0.158),
      env_double("V0_BEAM_Y", 0.285),
      env_double("V0_BEAM_Z", 0.0));

  const std::string fit_mode = env_string("V0_FIT_MODE", "kalman");
  if (!v0->set_track_fit_method(fit_mode))
  {
    std::cerr << "Unknown V0_FIT_MODE='" << fit_mode << "'" << std::endl;
    return 2;
  }

  v0->set_point_order(env_string("V0_POINT_ORDER", "radius"));
  v0->set_min_points(env_int("V0_MIN_POINTS", 20));
  v0->set_fit_first_points(env_int("V0_FIT_FIRST_POINTS", 8));
  v0->set_bfield(env_double("V0_BFIELD_T", 1.4));

  v0->use_kalman_field_map(env_bool("V0_USE_FIELD_MAP", false));
  v0->set_kalman_analytic_uniform_propagation(env_bool("V0_ANALYTIC_UNIFORM", true));
  v0->set_kalman_rkn4(
      env_double("V0_RKN_MAX_STEP_CM", 5.0),
      env_double("V0_RKN_TOLERANCE", 1.0e-4),
      env_int("V0_RKN_MAX_TRIALS", 12),
      env_int("V0_RKN_MAX_STEPS", 2000));

  // Deterministic, helix-like Kalman fit.
  v0->set_kalman_process_sigmas(
      env_double("V0_PROCESS_SIGMA_POS", 0.0),
      env_double("V0_PROCESS_SIGMA_PHI", 0.0),
      env_double("V0_PROCESS_SIGMA_QOPT", 0.0),
      env_double("V0_PROCESS_SIGMA_TANL", 0.0));
  v0->set_kalman_material(
      env_double("V0_MATERIAL_X0_PER_CM", 0.0),
      env_double("V0_SCATTERING_SCALE", 0.0),
      env_double("V0_ELOSS_GEV_PER_CM", 0.0),
      env_double("V0_ELOSS_SIGMA_FRACTION", 0.0));

  // Region-dependent measurement uncertainties [cm].
  v0->set_tpc_region_measurement_sigmas(
      env_double("V0_R1_SIGMA_RPHI", 0.040),
      env_double("V0_R1_SIGMA_R", 0.200),
      env_double("V0_R2_SIGMA_RPHI", 0.025),
      env_double("V0_R2_SIGMA_R", 0.050),
      env_double("V0_R3_SIGMA_RPHI", 0.025),
      env_double("V0_R3_SIGMA_R", 0.200),
      env_double("V0_SIGMA_Z", 0.100));
  v0->set_tpc_transition_measurement_sigmas(
      env_double("V0_TRANSITION_SIGMA_RPHI", 0.200),
      env_double("V0_TRANSITION_SIGMA_R", 0.500),
      env_double("V0_TRANSITION_SIGMA_Z", 0.200));
  v0->set_exclude_tpc_transition_layers(env_bool("V0_EXCLUDE_TRANSITION_LAYERS", true));

  // Fast common preselection. Keep this loose because it is evaluated before
  // the refined helix/Kalman pair PCA.
  v0->set_pre_track_pt_min(env_double("V0_PRE_TRACK_PT_MIN", 0.20));
  v0->set_pre_track_dca_xy_min(env_double("V0_PRE_TRACK_DCA_XY_MIN", -1.0));
  v0->set_pre_track_dca_z_min(env_double("V0_PRE_TRACK_DCA_Z_MIN", -1.0));
  v0->set_pre_track_dca_xy_max(env_double("V0_PRE_TRACK_DCA_XY_MAX", -1.0));
  v0->set_pre_track_dca_z_max(env_double("V0_PRE_TRACK_DCA_Z_MAX", -1.0));
  v0->set_pre_pair_dca_max(env_double("V0_PRE_PAIR_DCA_MAX", -1.0));
  v0->set_pre_lproj_min(env_double("V0_PRE_LPROJ_MIN", -1.0));
  v0->set_pre_cos_theta_min(env_double("V0_PRE_COSTHETA_MIN", -2.0));
  v0->set_pre_track_quality_max(env_double("V0_PRE_TRACK_QUALITY_MAX", -1.0));
  v0->set_pre_track_npoints_min(env_int("V0_PRE_TRACK_NPOINTS_MIN", 20));
  v0->set_pair_pca_z_max(env_double("V0_PAIR_PCA_Z_MAX", 20.0));

  // K0S -> pi+ pi-. These defaults reproduce the cuts you specified.
  v0->set_kshort_selection(
      env_bool("V0_ENABLE_KSHORT", true),
      env_int("V0_KSHORT_MIN_TPC_CLUSTERS", 20),
      env_double("V0_KSHORT_TRACK_PT_MIN", 0.20),
      env_double("V0_KSHORT_MASS_MIN", 0.40),
      env_double("V0_KSHORT_MASS_MAX", 0.60),
      env_double("V0_KSHORT_PCA_Z_MAX", 20.0),
      env_double("V0_KSHORT_PCA_DZ_MAX", 0.7),
      env_double("V0_KSHORT_DECAY_R_MIN", 2.0),
      env_double("V0_KSHORT_ALPHA_ABS_MAX", 1.0),
      env_double("V0_KSHORT_PAIR_DCA_MAX", 3.0),
      env_double("V0_KSHORT_DIRA_MIN", 0.75));

   v0->set_kshort_mass_range(
      env_double("V0_KSHORT_MASS_MIN", 0.40),
      env_double("V0_KSHORT_MASS_MAX", 0.60));

  // Lambda -> p pi- and anti-Lambda -> anti-p pi+.
  v0->set_lambda_selection(
      env_bool("V0_ENABLE_LAMBDA", true),
      env_int("V0_LAMBDA_MIN_TPC_CLUSTERS", 20),
      env_double("V0_LAMBDA_TRACK_PT_MIN", 0.20),
      env_double("V0_LAMBDA_MASS_MIN", 1.08),
      env_double("V0_LAMBDA_MASS_MAX", 1.16),
      env_double("V0_LAMBDA_PCA_Z_MAX", 20.0),
      env_double("V0_LAMBDA_PCA_DZ_MAX", 0.7),
      env_double("V0_LAMBDA_DECAY_R_MIN", 2.0),
      env_double("V0_LAMBDA_ALPHA_ABS_MAX", 1.0),
      env_double("V0_LAMBDA_PAIR_DCA_MAX", 3.0),
      env_double("V0_LAMBDA_DIRA_MIN", 0.75));

  v0->set_antilambda_selection(
      env_bool("V0_ENABLE_ANTILAMBDA", true),
      env_int("V0_ANTILAMBDA_MIN_TPC_CLUSTERS", 20),
      env_double("V0_ANTILAMBDA_TRACK_PT_MIN", 0.20),
      env_double("V0_ANTILAMBDA_MASS_MIN", 1.08),
      env_double("V0_ANTILAMBDA_MASS_MAX", 1.16),
      env_double("V0_ANTILAMBDA_PCA_Z_MAX", 20.0),
      env_double("V0_ANTILAMBDA_PCA_DZ_MAX", 0.7),
      env_double("V0_ANTILAMBDA_DECAY_R_MIN", 2.0),
      env_double("V0_ANTILAMBDA_ALPHA_ABS_MAX", 1.0),
      env_double("V0_ANTILAMBDA_PAIR_DCA_MAX", 3.0),
      env_double("V0_ANTILAMBDA_DIRA_MIN", 0.75));

  // Prompt phi -> K+K-. No secondary-vertex DIRA/decay-radius cut.
  v0->set_phi_selection(
    env_bool("V0_ENABLE_PHI", true),
    env_int("V0_PHI_MIN_TPC_CLUSTERS", 20),
    env_double("V0_PHI_TRACK_PT_MIN", 0.20),
    env_double("V0_PHI_MASS_MIN", 0.98),
    env_double("V0_PHI_MASS_MAX", 1.08),
    env_double("V0_PHI_TRACK_DCA_XY_MAX", 3.0),
    env_double("V0_PHI_TRACK_DCA_Z_MAX", -1.0),
    env_double("V0_PHI_PAIR_DCA_MAX", -1.0),
    env_double("V0_PHI_FLIGHT_LENGTH_MAX", -1.0),
    env_double("V0_PHI_PRIMARY_PCA_Z_MAX", 20.0),
    env_double("V0_PHI_PRIMARY_PCA_DZ_MAX", 1.0)
  );

  // Prompt D0 -> K- pi+ and anti-D0 -> K+ pi-.
  v0->set_d0_selection(
    env_bool("V0_ENABLE_D0", true),
    env_int("V0_D0_MIN_TPC_CLUSTERS", 20),
    env_double("V0_D0_TRACK_PT_MIN", 0.20),
    env_double("V0_D0_MASS_MIN", 1.70),
    env_double("V0_D0_MASS_MAX", 2.05),
    env_double("V0_D0_TRACK_DCA_XY_MAX", 3.0),
    env_double("V0_D0_TRACK_DCA_Z_MAX", -1.0),
    env_double("V0_D0_PAIR_DCA_MAX", -1.0),
    env_double("V0_D0_FLIGHT_LENGTH_MAX", -1.0),
    env_double("V0_D0_PRIMARY_PCA_Z_MAX", 20.0),
    env_double("V0_D0_PRIMARY_PCA_DZ_MAX", 2.0)
);

  v0->set_antid0_selection(
    env_bool("V0_ENABLE_ANTID0", true),
    env_int("V0_ANTID0_MIN_TPC_CLUSTERS", 20),
    env_double("V0_ANTID0_TRACK_PT_MIN", 0.20),
    env_double("V0_ANTID0_MASS_MIN", 1.70),
    env_double("V0_ANTID0_MASS_MAX", 2.05),
    env_double("V0_ANTID0_TRACK_DCA_XY_MAX", 3.0),
    env_double("V0_ANTID0_TRACK_DCA_Z_MAX", -1.0),
    env_double("V0_ANTID0_PAIR_DCA_MAX", -1.0),
    env_double("V0_ANTID0_FLIGHT_LENGTH_MAX", -1.0),
    env_double("V0_ANTID0_PRIMARY_PCA_Z_MAX", 20.0),
    env_double("V0_ANTID0_PRIMARY_PCA_DZ_MAX", 2.0));

  // Full daughter cluster payload remains K0S-only.
  v0->set_pair_pca_z_max(env_double("V0_PAIR_PCA_Z_MAX", 20.0));
  v0->set_write_same_sign_pairs(env_bool("V0_WRITE_SAME_SIGN", true));
  //v0->set_use_primary_vertex_kinematics(env_bool("V0_USE_PRIMARY_VERTEX_KINEMATICS", false));
  v0->set_write_track_tree(env_bool("V0_WRITE_TRACK_TREE", false));
  v0->set_write_cluster_residual_tree(env_bool("V0_WRITE_CLUSTER_RESIDUAL_TREE", false));
  v0->set_write_kalman_innovation_diagnostics(env_bool("V0_WRITE_KALMAN_DIAGNOSTICS", true));
  v0->set_write_kshort_daughter_details(env_bool("V0_WRITE_KSHORT_DETAILS", true));
  v0->set_kshort_detail_mass_range(
      env_double("V0_KSHORT_DETAIL_MASS_MIN", 0.40),
      env_double("V0_KSHORT_DETAIL_MASS_MAX", 0.60));
  v0->set_print_timing(env_bool("V0_PRINT_TIMING", false));
  v0->set_use_primary_vertex_kinematics(false);
  v0->set_write_same_sign_pairs(true);

  se->registerSubsystem(v0);

  auto *in = new Fun4AllDstInputManager("DSTin");
  if (in->AddListFile(input_list) != 0)
  {
    std::cerr << "Failed to open input list: " << input_list << std::endl;
    return 3;
  }
  se->registerInputManager(in);

  if (skipEvents > 0)
  {
    se->skip(skipEvents);
  }

  se->run(nEvents);
  se->End();
  delete se;

  std::cout << "Wrote " << output_file << std::endl;
  gSystem->Exit(0);
  return 0;
}

#endif

