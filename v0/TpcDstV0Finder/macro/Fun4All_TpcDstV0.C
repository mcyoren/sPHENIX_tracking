#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <TSystem.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libinmoduletracks.so)
R__LOAD_LIBRARY(libTpcDstV0Finder.so)

#include <TpcDstV0Finder.h>

void Fun4All_TpcDstV0(const char* input_dst,
                      const char* output_root = "TpcDstV0.root",
                      int nEvents = 0)
{
  auto* se = Fun4AllServer::instance();
  se->Verbosity(1);

  auto* v0 = new TpcDstV0Finder("TpcDstV0Finder", output_root);
  v0->set_track_node("FINALTRACKS");
  v0->set_vertex_node("FINALTRACKVERTICES");
  v0->use_vertex_node(true);          // falls back to (0,0,0) if node is absent
  v0->set_primary_vertex(0, 0, 0);
  v0->set_bfield(1.4);
  v0->write_track_tree(true);
  v0->write_same_sign_pairs(false);

  // Deliberately loose defaults for first validation.
  v0->set_min_clusters(20);
  v0->set_min_pt(0.15);
  v0->set_max_pair_dca(6.0);
  v0->set_min_decay_radius(0.0);
  v0->set_max_abs_decay_z(110.0);
  v0->set_min_lproj(0.0);
  v0->set_min_dira(-1.0);
  v0->set_max_abs_alpha(1.1);
  v0->set_search_range(120.0, 120.0);
  v0->set_coarse_steps(40);

  // Placeholder only. No refit is performed in this release.
  v0->enable_vertex_constrained_refit(false);
  se->registerSubsystem(v0);

  auto* in = new Fun4AllDstInputManager("DSTin");
  in->AddFile(input_dst);
  se->registerInputManager(in);

  se->run(nEvents);
  se->End();
  delete se;
  gSystem->Exit(0);
}
