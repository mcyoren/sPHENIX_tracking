#!/bin/bash
# space-charge map x keff scan (11 keff points per map) for one TPC geometry, local coordinates,
# one parallel job per map
# usage (from the macro directory): ./run_mapscan_geom.sh <g4geom> [outbase]
#   g4geom: tilt_0.16_0.85 (fitted tilt), survey, perfect, ...
MACRODIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GEOM=${1:?usage: run_mapscan_geom.sh <g4geom> [outbase]}
OUTBASE=${2:-/sphenix/tg/tg01/hf/mitrankov/cm_configs}
GTAG=$(echo ${GEOM} | tr . p)
OUT=${OUTBASE}/mapscan_${GTAG}_local
MAPS="cdb2d alpha1p25 alpha1p6 alpha2p0"
M=CMDistortionMaps_100ADC_
cd ${MACRODIR}

missing=0
for f in $(awk '!/^#/ && NF==8 {print $7; print $8}' mapscan/cm_mapscan_*_all.txt | grep -v '^-$' | sort -u); do
  [ -f "$f" ] || { echo "missing map: $f"; missing=1; }
done
[ $missing = 0 ] || exit 1

for m in ${MAPS}; do
  mkdir -p ${OUT}/${m}
  ./CMLaserAggregatedClustering_configs.sh 79516 0 ana571_2026p002_v001 run3pp ${OUT}/${m} ${GEOM} mapscan/cm_mapscan_${m}_all.txt 0 0 > ${OUT}/${m}/run.log 2>&1 &
  echo "started ${GEOM} ${m} (pid $!)"
done
wait

: > ${OUT}/all.list
for m in ${MAPS}; do
  echo "${m}: $(grep -c 'All done' ${OUT}/${m}/run.log) finished, $(ls ${OUT}/${m}/${M}*.root 2>/dev/null | wc -l) maps"
  awk '!/^#/ && NF==8 {print $1}' mapscan/cm_mapscan_${m}_all.txt | while read tag; do
    f=${OUT}/${m}/${M}${tag}_g4${GTAG}_local.root
    [ -f "$f" ] && echo "$f ${tag}" >> ${OUT}/all.list
  done
done
grep -v noGarfield ${OUT}/all.list > ${OUT}/keff_fit.list
root -l -b -q "FitCMKeff.C+(\"${OUT}/keff_fit.list\",\"${OUT}/keff_fit.pdf\")" 2>&1 | grep -v "^Info\|^$\|Processing" | tee ${OUT}/keff_fit.txt
