#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_tpc_dv_scan.sh runlist.txt /path/to/Fun4All_TpcDriftVelocityScan.C output_dir
#
# runlist.txt columns:
#   run_number   one_DST_file_for_that_run   [measured_t0_in_units_of_8_timebins]
#
# Example:
#   79513 /path/to/run79513_one_segment.root  -999
#   79522 /path/to/run79522_one_segment.root  0.125
#
# One independent ROOT/Fun4All process is launched per run. This avoids carrying
# Fun4AllServer/CDB singleton state from one run into the next.

RUNLIST=${1:?need run list}
MACRO=${2:-Fun4All_TpcDriftVelocityScan.C}
OUTDIR=${3:-dvscan}
NPHI=${NPHI:-24}
STEP_NS=${STEP_NS:-56.8}
GLOBALTAG=${GLOBALTAG:-newcdbtag}

mkdir -p "$OUTDIR"

while read -r run input t0 rest; do
  [[ -z "${run:-}" ]] && continue
  [[ "$run" =~ ^# ]] && continue
  t0=${t0:--999}

  echo "=== run $run ==="
  root -l -b -q \
    "${MACRO}(${run},\"${input}\",${t0},${NPHI},${STEP_NS},\"${OUTDIR}\",\"${GLOBALTAG}\")"
done < "$RUNLIST"
