#ifndef MACRO_FUN4ALLLASERAGGREGATEDCLUSTERINGCONFIGS_C
#define MACRO_FUN4ALLLASERAGGREGATEDCLUSTERINGCONFIGS_C

// Ben's Fun4All_LaserAggregatedClustering_simple.C, extended so that one job
// clusters once and then runs one LaserClusterTruthMatcher per PHGarfield
// configuration listed in a config file. The G4/Acts TPC geometry
// (survey / perfect / as stored in the CDB geometry file) is one per job.
//
// config file: one matcher per line, '#' = comment, '-' = default
//   tag  garfield(0|1)  tpcgeom(hard|container|perfect)  keff0  keff1  efield2D  efield3D_side0  efield3D_side1
// e.g.
//   noGarfield   0  -          -    -    -  -  -
//   hard_cdb     1  hard       -    -    -  -  -
//   cont_keff1   1  container  1.0  1.0  -  -  -
//   cont_k105    1  container  x1.05 x1.05 -  -  -   ('x' = factor times the CDB keff)

#include <GlobalVariables.C>
#include <G4_ActsGeom.C>
#include <Trkr_RecoInit.C>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <phool/recoConsts.h>

#include <tpc/LaserAggregatedClusterizer.h>
#include <tpc/LaserClusterHelper.h>

#include <tpccalib/LaserClusterTruthMatcher.h>
#include <tpccalib/TpcGeometryOverride.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libtpccalib.so)

namespace CMCONFIG
{
  struct MatcherConfig
  {
    std::string tag;
    bool garfield{true};
    std::string tpcgeom{"hard"};
    double keff0{-1};  // <0: from CDB
    double keff1{-1};
    double keffScale0{-1};  // >0: factor times the CDB keff
    double keffScale1{-1};
    std::string efield2D;
    std::string efield3D0;
    std::string efield3D1;
  };

  std::vector<MatcherConfig> readConfigs(const std::string &configfile)
  {
    std::vector<MatcherConfig> configs;
    std::ifstream in(configfile);
    std::string line;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      std::istringstream ss(line);
      std::string tag, garf, geom, k0, k1, e2, e30, e31;
      if (!(ss >> tag >> garf >> geom >> k0 >> k1 >> e2 >> e30 >> e31))
      {
        std::cout << "CMCONFIG - skipping malformed line: " << line << std::endl;
        continue;
      }
      MatcherConfig c;
      c.tag = tag;
      c.garfield = (garf == "1");
      if (geom != "-") { c.tpcgeom = geom; }
      if (k0 != "-") { (k0[0] == 'x' ? c.keffScale0 : c.keff0) = std::stod(k0[0] == 'x' ? k0.substr(1) : k0); }
      if (k1 != "-") { (k1[0] == 'x' ? c.keffScale1 : c.keff1) = std::stod(k1[0] == 'x' ? k1.substr(1) : k1); }
      if (e2 != "-") { c.efield2D = e2; }
      if (e30 != "-") { c.efield3D0 = e30; }
      if (e31 != "-") { c.efield3D1 = e31; }
      configs.push_back(c);
    }
    return configs;
  }
}  // namespace CMCONFIG

//! g4geom: "cdb" = TPC place/rot as stored in Tracking_Geometry, "survey" or "perfect" = overridden,
//!         "tilt_<rot_x>_<rot_y>" = survey translation with custom tilt in mrad
void Fun4All_LaserAggregatedClustering_configs(const std::string &listfile = "segments.list",
                                               const std::string &outdir = ".",
                                               const std::string &outfile = "LaserAggregatedClusters.root",
                                               const std::string &g4geom = "cdb",
                                               const std::string &configfile = "cm_garfield_configs.txt",
                                               const bool doQA = false,
                                               const bool useGlobal = false)  // local TPC frame: truth pattern is in TPC coordinates
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  // ---- pull the run number off the first entry in the list -----------------
  int runnumber = 0;
  {
    std::ifstream in(listfile);
    std::string firstline;
    if (!in.is_open() || !std::getline(in, firstline))
    {
      std::cout << "could not read " << listfile << std::endl;
      gSystem->Exit(1);
      return;
    }
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstline);
    runnumber = runseg.first;
  }
  std::cout << "runnumber: " << runnumber << std::endl;

  const auto configs = CMCONFIG::readConfigs(configfile);
  if (configs.empty())
  {
    std::cout << "no matcher configs in " << configfile << std::endl;
    gSystem->Exit(1);
    return;
  }

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  CDBInterface::instance()->Verbosity(1);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  // ---- G4/Acts TPC geometry: must come before TrackingInit (MakeActsGeometry)
  if (g4geom != "cdb")
  {
    TpcGeometryOverride *tpcgeo = new TpcGeometryOverride;
    if (g4geom == "perfect")
    {
      tpcgeo->set_perfect();
    }
    else if (g4geom == "survey")
    {
      tpcgeo->set_survey();
    }
    else if (g4geom.rfind("tilt_", 0) == 0)
    {
      // tilt_<rot_x>_<rot_y> in mrad, survey translation, e.g. tilt_-0.15_0.74
      double rx = 0;
      double ry = 0;
      if (std::sscanf(g4geom.c_str(), "tilt_%lf_%lf", &rx, &ry) != 2)
      {
        std::cout << "cannot parse " << g4geom << " (tilt_<rot_x mrad>_<rot_y mrad>)" << std::endl;
        gSystem->Exit(1);
        return;
      }
      tpcgeo->set_survey();
      tpcgeo->set_rot(rx * 1e-3, ry * 1e-3, 0);
    }
    else
    {
      std::cout << "unknown g4geom " << g4geom << " (cdb|survey|perfect|tilt_<rx>_<ry>)" << std::endl;
      gSystem->Exit(1);
      return;
    }
    se->registerSubsystem(tpcgeo);
  }

  TRACKING::streaming_mode = true;

  TrackingInit();

  // ---- the clusterizer: merges the run trees in InitRun, clusters once -----
  LaserAggregatedClusterizer *clus = new LaserAggregatedClusterizer("LaserAggregatedClusterizer");
  clus->Verbosity(0);
  if (doQA)
  {
    clus->set_QAFile(outdir + "/" + "LASER_AGGREGATED_CLUSTER_100ADC_QA.root");
  }
  clus->set_padContainerNodeName("LASER_AGGREGATED_PADS_100ADC");
  clus->set_clusterNodeName("LASER_AGGREGATED_CLUSTER_100ADC");
  clus->SetSegmentList(listfile);  // must precede registerSubsystem
  se->registerSubsystem(clus);

  // ---- CDB keff, for configs that scale it --------------------------------
  double cdbKeff[2] = {1.0, 1.0};
  {
    const std::string kefffile = CDBInterface::instance()->getUrl("Tpc_PolyClusterizer_kEff");
    if (!kefffile.empty())
    {
      CDBTTree keffTree(kefffile);
      keffTree.LoadCalibrations();
      cdbKeff[0] = keffTree.GetSingleFloatValue("keffside0");
      cdbKeff[1] = keffTree.GetSingleFloatValue("keffside1");
    }
    std::cout << "CDB keff: side0 = " << cdbKeff[0] << ", side1 = " << cdbKeff[1] << std::endl;
  }

  // ---- one matcher per config ---------------------------------------------
  for (const auto &c : configs)
  {
    std::string name = c.tag + "_g4" + g4geom + (useGlobal ? "" : "_local");
    std::replace(name.begin(), name.end(), '.', 'p');  // Fun4All module names may not contain dots

    LaserClusterTruthMatcher *matcher = new LaserClusterTruthMatcher("LaserClusterTruthMatcher_" + name);
    matcher->Verbosity(0);
    matcher->setLaserClusterNodeName("LASER_AGGREGATED_CLUSTER_100ADC");
    matcher->setLaserClusterNodeNameOut("LASER_AGGREGATED_CLUSTER_100ADC_MATCHED_" + name);
    matcher->setTruthFile("/sphenix/u/bkimelman/CMStripePattern_full_Aug19_2026.root");
    matcher->set_useGlobal(useGlobal);
    matcher->set_usePHGarfieldDistortions(c.garfield);
    if (c.garfield)
    {
      if (c.tpcgeom == "hard")
      {
        matcher->set_garfield_tpc_geometry(LaserClusterHelper::GarfieldTpcGeometry::HardcodedSurvey);
      }
      else if (c.tpcgeom == "container")
      {
        matcher->set_garfield_tpc_geometry(LaserClusterHelper::GarfieldTpcGeometry::GeomContainer);
      }
      else if (c.tpcgeom == "perfect")
      {
        matcher->set_garfield_tpc_geometry(LaserClusterHelper::GarfieldTpcGeometry::Perfect);
      }
      else
      {
        std::cout << "unknown tpcgeom " << c.tpcgeom << " for " << c.tag << std::endl;
        gSystem->Exit(1);
        return;
      }
      if (c.keff0 >= 0) { matcher->set_garfield_keffside0(c.keff0); }
      if (c.keff1 >= 0) { matcher->set_garfield_keffside1(c.keff1); }
      if (c.keffScale0 > 0) { matcher->set_garfield_keffside0(c.keffScale0 * cdbKeff[0]); }
      if (c.keffScale1 > 0) { matcher->set_garfield_keffside1(c.keffScale1 * cdbKeff[1]); }
      if (!c.efield2D.empty()) { matcher->set_garfield_efieldmap(c.efield2D); }
      if (!c.efield3D0.empty() && !c.efield3D1.empty()) { matcher->set_garfield_efieldmap3d(c.efield3D0, c.efield3D1); }
    }
    if (doQA)
    {
      matcher->set_QABase(outdir + "/" + "QA_LASER_AGGREGATED_Matching_100ADC_" + name);
    }
    matcher->setDistortionOutputfile(outdir + "/" + "CMDistortionMaps_100ADC_" + name + ".root");
    matcher->set_grid_dimensions(36, 16);
    se->registerSubsystem(matcher);
  }

  // ---- no event data is read; the dummy manager just drives one event ------
  Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(in);

  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outdir + "/" + outfile);
  se->registerOutputManager(out);

  se->run(1);
  se->End();

  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
}

#endif
