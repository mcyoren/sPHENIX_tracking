DC synchronized preprocessing, first version
============================================

Put these files in:
  TrackingAnalysis/DC_ML_processing/1_dc_preprocessing/

Keep the existing files there too:
  T_DigitalCurrent.C
  T_DigitalCurrent.h
  InfoTPC.h
  R1_ChannelMapping.csv
  R2_ChannelMapping.csv
  R3_ChannelMapping.csv

Files in this package:
  build_dc_index.C
  make_dc_windows.C
  make_dc_maps.C
  run_dc_window.sh
  dc_windows.sub
  prepare_and_submit_dc.sh

Quick run for 81558 with 10 synchronized windows
-------------------------------------------------

  ./prepare_and_submit_dc.sh 81558 10 \
    /sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run81558

The script will:
  1) index the BCO coverage of every DC ROOT file;
  2) find the common coverage of all 24 EBDCs x 2 streams;
  3) divide the common live coverage into 10 equal-live-frame windows;
  4) submit one Condor job per window.

Each output window keeps four versions:
  raw             = dc_current
  ped60_signed    = dc_current - 60*dc_nsamples, negatives kept
  ped60_positive  = max(0, dc_current - 60*dc_nsamples)
  ped60_strict    = existing preprocessing choice: 0 <= ped60 <= 30000

No gain correction is applied here.

The output also contains:
  - side-separated phi x r x frame TH3 maps for every version;
  - side-separated phi x r x z TH3 maps for every version;
  - native pad x layer maps for side/sector/R1,R2,R3 for every version;
  - nsamples exposure maps;
  - sample-block count maps;
  - h_dc_norm with additive normalization quantities;
  - h_dc_rejection_stats;
  - dc_metadata tree.

Notebook-compatible aliases currently point to ped60_strict:
  h_adc_pad_vs_layer_side{side}_sec{sector:02d}_R{module}_PHGarfieldRawHitsQA
  h_adc_phi_vs_radius_vs_tbin_side{side}_PHGarfieldRawHitsQA
  h_nEvents_PHGarfieldRawHitsQA

For DC input, the notebook time axis is 80 DC frame bins, not the old ~970 raw-hit
TPC time bins. For the first check, use the complete DC frame range rather than the old
trigger-crossing time cuts.

After jobs finish
-----------------

  hadd -f dc_run81558/dc_run81558_all.root \
    dc_run81558/output/dc_run00081558_win*.root

h_dc_norm is intentionally additive under hadd. In particular it stores, per side:
  N_WINDOWS
  LIVE_BCO_TICKS
  LIVE_FRAME_SLOTS
  LIVE_TPC_VOLUMES
  CRC_GOOD_ENTRIES
  MAPPED_SAMPLE_BLOCKS
  NSAMPLES_SUM
  RAW_CURRENT_SUM
  PED60_SIGNED_SUM
  PED60_POSITIVE_SUM
  PED60_STRICT_SUM
  OBSERVED_UNIQUE_FRAME_SLOTS

This lets us later define and compare DC normalization observables for kEff without
rerunning the trees.

Important first-pass caveats
----------------------------

1) The fixed pedestal 60 is copied from the existing T_DigitalCurrent.C. We keep other
   map versions so we can change the choice later without losing the raw information.
2) There is not enough information in the current DC tree to identify every individual
   saturated SAMPA ADC sample. The strict >30000 rejection is retained as a QA/cleaning
   variant, not claimed as a complete saturation correction.
3) The z/frame phase uses floor(gtm_bco/9720) mod 80. The histogram contents are useful
   now, but the absolute z phase should be validated before using the 3D map as a final
   PHGarfield source.

============================================================
R1 GEM LOAD CURRENT FROM TPC_CONDITIONS
============================================================

Each output window ROOT file now also queries TPC_CONDITIONS using the run
number and one representative BCO at the midpoint of the window's accumulated
live exposure.

The standalone helper DCGemCurrent.h reproduces the SR1/NR1 channel lists and
the interpolated-median-current logic used by TpcConditionsReco.

Stored per window in dc_metadata:
  conditions_available
  conditions_bco
  conditions_nearest_bco
  conditions_first_bco
  conditions_last_bco
  gem_load_SR1
  gem_load_NR1
  gem_average_SR1
  gem_average_NR1
  cdb_global_tag
  tpc_conditions_url

Side convention:
  side0 = South = LoadSR1
  side1 = North = LoadNR1

h_dc_norm has four additional hadd-safe quantities:
  GEM_R1_CURRENT_SUM
  GEM_R1_CURRENT2_SUM
  GEM_R1_CURRENT_X_LIVE_FRAMES
  GEM_R1_CURRENT_SAMPLES

For a hadded run file:
  unweighted current = GEM_R1_CURRENT_SUM / GEM_R1_CURRENT_SAMPLES
  live-frame-weighted current = GEM_R1_CURRENT_X_LIVE_FRAMES / LIVE_FRAME_SLOTS

p_dc_gem_r1_current is a TProfile and directly gives the live-frame-weighted
R1 GEM current after hadd.

The default CDB tag is newcdbtag, matching the current TPC production macro.
Override it with the fourth argument to prepare_and_submit_dc.sh if needed.

Quick login-node test after windows are made:
  root -l -b -q 'test_dc_gem_current.C(81558,957442302614,"newcdbtag")'

Submit all Takao runs with 10 windows/run:
  ./submit_all_takao_runs.sh 10 newcdbtag
