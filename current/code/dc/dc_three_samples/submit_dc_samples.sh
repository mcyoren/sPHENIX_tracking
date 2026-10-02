#!/usr/bin/env bash
set -euo pipefail

# Submit the already-prepared 3 sample windows for one run.
#
# Usage:
#   ./submit_dc_samples.sh RUN [CDB_GLOBALTAG]

RUN="${1:?run number required}"
CDBTAG="${2:-newcdbtag}"

SRCDIR="$(cd "$(dirname "$0")" && pwd)"
WORKDIR="$SRCDIR/dc_samples_run${RUN}"
OUTDIR="$WORKDIR/output"
mkdir -p "$OUTDIR"

OLDINDEX="$SRCDIR/dc_run${RUN}/dc_index_run${RUN}.csv"
INDEX="$WORKDIR/dc_index_run${RUN}.csv"
if [ -s "$OLDINDEX" ]; then
  INDEX="$OLDINDEX"
fi

WINDOWS="$WORKDIR/dc_sample_windows_run${RUN}.csv"
JOBS="$WORKDIR/dc_sample_windows_run${RUN}.jobs"

for f in "$INDEX" "$WINDOWS" "$JOBS"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing $f"
    echo "Run prepare_dc_samples.sh first."
    exit 2
  fi
done

echo "Submitting 3 samples for run $RUN"
condor_submit "$SRCDIR/dc_windows.sub" \
  run="$RUN" \
  index="$INDEX" \
  windows="$WINDOWS" \
  jobs="$JOBS" \
  outdir="$OUTDIR" \
  srcdir="$SRCDIR" \
  mappingdir="$SRCDIR" \
  cdbtag="$CDBTAG"

printf -v RUN8 "%08d" "$RUN"
echo
echo "After the 3 jobs finish:"
echo "  hadd -f ${WORKDIR}/dc_run${RUN}_3samples.root ${OUTDIR}/dc_run${RUN8}_win*.root"
