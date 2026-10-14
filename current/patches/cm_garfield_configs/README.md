# CM laser distortion maps: PHGarfield configurations, TPC geometry, keff and map scans

On top of bkimelman:laserCluster_hitChange (coresoftware PR 4480).

## coresoftware
- `0001-coresoftware-cm-garfield-configs.patch` (git am): PHGarfield field map (2D / 3D), keff and
  TPC placement (hard / container / perfect) settable from LaserClusterTruthMatcher; private
  PHGarfield per matcher; `TpcGeometryOverride` module (TPC place/rot in TPCGEOMCONTAINER).

## Running (macros/TrackingProduction)
- `Fun4All_LaserAggregatedClustering_configs.C`: clusters once, one matcher per config line.
  Arguments: list, outdir, outfile, g4geom (`cdb|survey|perfect|tilt_<rx>_<ry>` in mrad), config
  file, doQA, useGlobal (default false: local TPC frame, the truth pattern is in TPC coordinates).
- `CMLaserAggregatedClustering_configs.sh`: job script for the macro.
- Config files: one matcher per line, `tag garfield tpcgeom keff0 keff1 efield2D efield3D_side0 efield3D_side1`,
  `-` = default, keff `x<f>` = f times the CDB keff.
  `cm_garfield_configs.txt`, `cm_garfield_keffscan.txt`, `cm_garfield_tiltscan.txt`, `mapscan/cm_mapscan_<map>[_ext|_all].txt`.
- `run_mapscan.sh`: map x keff scan, ideal TPC (global coordinates, as first done).
- `run_mapscan_geom.sh <g4geom>`: map x keff scan (11 keff points) for one TPC geometry in local
  coordinates, e.g. `tilt_0.16_0.85` (fitted) or `survey`.

## Analysis
- `CMQA.C`: multi-page QA of a list of maps (Ben-style projections, 2D maps, differences,
  higher-order residuals, vs phi / vs R, per-ring decomposition, summary).
- `FitCMKeff.C`: keff fit per map and side from a keff scan (ring means, M0 / M1 with nuisances).
- `PlotCMConfigProjections.C`, `PlotCMSurveyCheck.C`: projections and scores; first-harmonic (tilt) fit.
- `make_qa.sh [base] [outdir]`: all QA studies from the outputs under `cm_configs`.
- `qa_run79516/`: numbers from the QA of run 79516 (pdfs not committed).

## Findings, run 79516
- Global coordinates vs a TPC-frame truth pattern turn the survey shift into a fake distortion
  (Ben: use local coordinates or transform the truth) -> `ben_message.md`.
- The laser sees the TPC tilt relative to B (electrons follow B), not a rigid shift: fitted
  rot_x = 0.16, rot_y = 0.85 mrad (survey -0.30, 1.48). Without tilt, North shows a ~1.3 mm first harmonic.
- keff fit with nuisances: South and North agree; alpha1p6 best at ~1.12-1.19 x CDB keff.
- Remaining, configuration independent: North ring-mean offsets (~+0.7 mm dR, ~1.2-1.7 mrad rotation)
  and a ring-to-ring dR zigzag of +-0.5-0.8 mm on both sides.
