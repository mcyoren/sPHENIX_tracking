#!/usr/bin/env bash
set -euo pipefail

NWINDOWS="${1:-10}"
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
  echo "======================================"
  echo "Submitting run $run with $NWINDOWS windows"
  echo "Input: $d"
  echo "CDB tag: $CDBTAG"
  echo "======================================"

  "$SRCDIR/prepare_and_submit_dc.sh" "$run" "$NWINDOWS" "$d" "$CDBTAG"
done
