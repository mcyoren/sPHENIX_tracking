#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./prepare_and_submit_dc.sh RUN NWINDOWS INPUTDIR [CDB_GLOBALTAG]
# Example:
#   ./prepare_and_submit_dc.sh 81558 10 /sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run81558 newcdbtag

RUN="${1:?run number required}"
NWINDOWS="${2:-10}"
INPUTDIR="${3:-/sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run${RUN}}"
CDBTAG="${4:-newcdbtag}"

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
MAPPINGDIR="$SRCDIR"
WORKDIR="$SRCDIR/dc_run${RUN}"
OUTDIR="$WORKDIR/output"
LOGDIR="$WORKDIR/logs"

mkdir -p "$WORKDIR" "$OUTDIR" "$LOGDIR"

INDEX="$WORKDIR/dc_index_run${RUN}.csv"
WINDOWS="$WORKDIR/dc_windows_run${RUN}.csv"
SUMMARY="$WORKDIR/dc_windows_summary_run${RUN}.csv"
JOBS="$WORKDIR/dc_windows_run${RUN}.jobs"

cd "$SRCDIR"
printf -v RUN8 "%08d" "$RUN"

echo "=== Building DC file/BCO index ==="
root -l -b -q \
  "build_dc_index.C(${RUN},\"${INPUTDIR}\",\"${INDEX}\")"

echo "=== Making ${NWINDOWS} synchronized BCO windows ==="
root -l -b -q \
  "make_dc_windows.C(${RUN},${NWINDOWS},\"${INDEX}\",\"${WINDOWS}\",\"${SUMMARY}\",\"${JOBS}\")"

echo "=== Window summary ==="
cat "$SUMMARY"

echo "=== Submitting Condor jobs ==="
echo "CDB tag: $CDBTAG"
condor_submit "$SRCDIR/dc_windows.sub" \
  run="$RUN" \
  index="$INDEX" \
  windows="$WINDOWS" \
  jobs="$JOBS" \
  outdir="$OUTDIR" \
  srcdir="$SRCDIR" \
  mappingdir="$MAPPINGDIR" \
  logdir="$LOGDIR" \
  cdbtag="$CDBTAG"

echo
echo "After jobs finish:"
echo "  hadd -f ${WORKDIR}/dc_run${RUN}_all.root ${OUTDIR}/dc_run${RUN8}_win*.root"
echo "DC normalization and GEM-current sums are preserved in h_dc_norm."
echo "Per-window GEM currents are preserved in dc_metadata."
echo "p_dc_gem_r1_current becomes the live-frame-weighted current after hadd."
