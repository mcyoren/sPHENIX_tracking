THREE SHORT DC SAMPLES PER RUN
==============================

Purpose
-------
Instead of processing a large fraction of every run, use three short,
synchronized samples:
  window 0: early  = 10% of common live exposure
  window 1: middle = 50%
  window 2: late   = 90%

The default production width is 2000 DC frame slots per sample, i.e.
25 complete 80-frame TPC volumes. This is intentionally small compared
with a full run but still gives a large number of DC channel measurements.

Files
-----
make_dc_sample_windows.C
prepare_dc_samples.sh
submit_dc_samples.sh
submit_all_takao_samples.sh

These reuse your existing:
  build_dc_index.C
  make_dc_maps.C
  run_dc_window.sh
  dc_windows.sub

Quick test on run 81558
-----------------------
1) Copy the four new files into:
   TrackingAnalysis/DC_ML_processing/1_dc_preprocessing/

2) Make scripts executable:
   chmod +x prepare_dc_samples.sh submit_dc_samples.sh submit_all_takao_samples.sh

3) Generate very small test samples (400 frames = 5 TPC volumes each):
   ./prepare_dc_samples.sh 81558 400 \
     /sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run81558 \
     newcdbtag

4) Run only the middle sample interactively:
   SRCDIR=$PWD
   mkdir -p dc_samples_run81558/test_output

   ./run_dc_window.sh \
     81558 \
     1 \
     "$SRCDIR/dc_run81558/dc_index_run81558.csv" \
     "$SRCDIR/dc_samples_run81558/dc_sample_windows_run81558.csv" \
     "$SRCDIR/dc_samples_run81558/test_output" \
     "$SRCDIR" \
     "$SRCDIR" \
     newcdbtag \
     2>&1 | tee /tmp/dc_81558_middle_small.log

If no full-run index exists, prepare_dc_samples.sh creates one in
dc_samples_run81558 and prints the exact interactive command using that path.

Production: all Takao runs
--------------------------
After the test looks good:

   ./submit_all_takao_samples.sh 2000 newcdbtag

This submits exactly 3 Condor jobs per run.

Merge the 3 outputs of one run
------------------------------
For run 81558:

   hadd -f dc_samples_run81558/dc_run81558_3samples.root \
     dc_samples_run81558/output/dc_run00081558_win*.root

h_dc_norm remains additive, and dc_metadata keeps one entry per sample.
