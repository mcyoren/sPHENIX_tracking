#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 8 ]; then
  echo "Usage: $0 RUN WINDOW INDEX_CSV WINDOWS_CSV OUTDIR SRCDIR MAPPINGDIR CDB_GLOBALTAG"
  exit 2
fi

RUN="$1"
WINDOW="$2"
INDEX_CSV="$3"
WINDOWS_CSV="$4"
OUTDIR="$5"
SRCDIR="$6"
MAPPINGDIR="$7"
CDBTAG="$8"

mkdir -p "$OUTDIR"
OUTFILE=$(printf "%s/dc_run%08d_win%02d.root" "$OUTDIR" "$RUN" "$WINDOW")

cd "$SRCDIR"

echo "Host: $(hostname)"
echo "Run: $RUN  window: $WINDOW"
echo "CDB global tag: $CDBTAG"
echo "Output: $OUTFILE"

time root -l -b -q \
  "make_dc_maps.C(${RUN},${WINDOW},\"${INDEX_CSV}\",\"${WINDOWS_CSV}\",\"${OUTFILE}\",\"${MAPPINGDIR}\",\"${CDBTAG}\")"
