#!/usr/bin/env bash
set -euo pipefail

# Prepare and submit 3 short synchronized samples for every Takao run.
#
# Defaults:
#   2000 frames/sample = 25 TPC volumes/sample
#   positions = 10%, 50%, 90% of common live exposure
#
# Usage:
#   ./submit_all_takao_samples.sh [SAMPLE_FRAMES] [CDB_GLOBALTAG] [BASE]
#
# Example:
#   ./submit_all_takao_samples.sh 2000 newcdbtag

SAMPLE_FRAMES="${1:-2000}"
CDBTAG="${2:-newcdbtag}"
BASE="${3:-/sphenix/lustre01/sphnxpro/t_sakagu/DCTrees}"
SRCDIR="$(cd "$(dirname "$0")" && pwd)"

for d in "$BASE"/run*; do
  [ -d "$d" ] || continue
  run=$(basename "$d" | sed 's/^run//')

  if ! compgen -G "$d/DCTree_DC-ebdc*.root" > /dev/null; then
    echo "Skipping run $run: no DC trees found"
    continue
  fi

  echo
  echo "=================================================="
  echo "Preparing run $run"
  echo "3 samples x $SAMPLE_FRAMES frames"
  echo "Input: $d"
  echo "=================================================="

  "$SRCDIR/prepare_dc_samples.sh" \
    "$run" "$SAMPLE_FRAMES" "$d" "$CDBTAG"

  "$SRCDIR/submit_dc_samples.sh" \
    "$run" "$CDBTAG"
done
