#!/usr/bin/env bash
set -euo pipefail

# Build/reuse the run index and generate 3 short synchronized samples:
#   window 0 = early  (10% of common live exposure)
#   window 1 = middle (50%)
#   window 2 = late   (90%)
#
# Usage:
#   ./prepare_dc_samples.sh RUN SAMPLE_FRAMES INPUTDIR [CDB_GLOBALTAG]
#
# Example quick test:
#   ./prepare_dc_samples.sh 81558 400 \
#     /sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run81558 newcdbtag
#
# Example production sampling:
#   ./prepare_dc_samples.sh 81558 2000 \
#     /sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run81558 newcdbtag

RUN="${1:?run number required}"
SAMPLE_FRAMES="${2:-2000}"
INPUTDIR="${3:-/sphenix/lustre01/sphnxpro/t_sakagu/DCTrees/run${RUN}}"
CDBTAG="${4:-newcdbtag}"

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
WORKDIR="$SRCDIR/dc_samples_run${RUN}"
OUTDIR="$WORKDIR/output"
mkdir -p "$WORKDIR" "$OUTDIR"

# Reuse an already-built full-run index if available.
OLDINDEX="$SRCDIR/dc_run${RUN}/dc_index_run${RUN}.csv"
INDEX="$WORKDIR/dc_index_run${RUN}.csv"

if [ -s "$OLDINDEX" ]; then
  INDEX="$OLDINDEX"
  echo "=== Reusing existing index ==="
  echo "$INDEX"
elif [ -s "$INDEX" ]; then
  echo "=== Reusing existing sample index ==="
  echo "$INDEX"
else
  echo "=== Building DC file/BCO index ==="
  root -l -b -q \
    "build_dc_index.C(${RUN},\"${INPUTDIR}\",\"${INDEX}\")"
fi

WINDOWS="$WORKDIR/dc_sample_windows_run${RUN}.csv"
SUMMARY="$WORKDIR/dc_sample_windows_summary_run${RUN}.csv"
JOBS="$WORKDIR/dc_sample_windows_run${RUN}.jobs"

echo "=== Making 3 short synchronized samples, ${SAMPLE_FRAMES} frames each ==="
root -l -b -q \
  "make_dc_sample_windows.C(${RUN},${SAMPLE_FRAMES},\"${INDEX}\",\"${WINDOWS}\",\"${SUMMARY}\",\"${JOBS}\")"

echo
echo "=== Sample summary ==="
cat "$SUMMARY"

echo
echo "Prepared only. Nothing submitted yet."
echo
echo "Quick manual test of the MIDDLE sample (window 1):"
echo
echo "  ./run_dc_window.sh \\"
echo "    ${RUN} 1 \\"
echo "    \"${INDEX}\" \\"
echo "    \"${WINDOWS}\" \\"
echo "    \"${WORKDIR}/test_output\" \\"
echo "    \"${SRCDIR}\" \\"
echo "    \"${SRCDIR}\" \\"
echo "    \"${CDBTAG}\""
