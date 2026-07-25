# TpcV0CandidateTree local and Condor runner

Files:

- `Fun4All_TpcV0CandidateTree.C`: configures `TpcV0CandidateTree` from environment variables.
- `run_tpc_v0_job.sh`: worker used locally and by HTCondor.
- `submit_tpc_v0.sh`: discovers DST files, creates file lists and either runs one local job or submits one Condor job per group of files.

The default fit is Kalman with zero process noise, zero material, and the requested R1/R2/R3 uncertainties. K0S daughter details are enabled only for the configurable K0S detail mass window; the all-track tree is disabled.

## Local test

```bash
./submit_tpc_v0.sh \
  --mode local \
  --input-dir /path/to/DST \
  --pattern 'HITS_pp_79513_*.root' \
  --output-dir ./output \
  --output-base v0_test \
  --nevents 100 \
  --pair-pca-z-max 10 \
  --pair-pca-dz-max 0.5 \
  --pair-decay-r-min 6 \
  --pair-dca-max 1.0 \
  --pair-dira-min 0.99 \
  --kshort-mass-min 0.44 \
  --kshort-mass-max 0.56
```

Local mode intentionally runs only the first generated list. Increase `--files-per-job` to place several input files in that one list.

## Condor production

```bash
./submit_tpc_v0.sh \
  --mode condor \
  --input-dir /path/to/DST \
  --pattern 'HITS_pp_79513_*.root' \
  --output-dir /path/to/output \
  --output-base v0_79513 \
  --files-per-job 2 \
  --memory 6GB \
  --pair-pca-z-max 10 \
  --pair-pca-dz-max 0.5 \
  --pair-decay-r-min 6 \
  --pair-dca-max 1.0 \
  --pair-dira-min 0.99
```

The generated work area is:

```text
output/condor_<output-base>/
  lists/
  log/
  v0_config.sh
  submit_<output-base>.sub
```

## Environment

Run from an already initialized sPHENIX environment, or pass a setup script:

```bash
--setup-script /absolute/path/to/setup_sphenix.sh
```

The setup script must make ROOT, Fun4All, and the installed `libTrackingDiagnostics.so` available.
