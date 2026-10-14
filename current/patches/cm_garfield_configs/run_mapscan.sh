#!/bin/bash
# space-charge map x keff scan, perfect TPC geometry, one parallel job per map
# usage (from the macro directory): ./run_mapscan.sh [outbase] [config suffix, e.g. _ext]
MACRODIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${1:-/sphenix/tg/tg01/hf/mitrankov/cm_configs/mapscan}
MAPS="cdb2d alpha1p25 alpha1p6 alpha2p0"
SUFFIX=${2:-}
cd ${MACRODIR}

# all 3D map files must exist before anything starts
missing=0
for f in $(awk '!/^#/ && NF==8 {print $7; print $8}' mapscan/cm_mapscan_*${SUFFIX}.txt | grep -v '^-$' | sort -u); do
  [ -f "$f" ] || { echo "missing map: $f"; missing=1; }
done
[ $missing = 0 ] || exit 1

for m in ${MAPS}; do
  mkdir -p ${OUT}/${m} ${OUT}/${m}${SUFFIX}
  ./CMLaserAggregatedClustering_configs.sh 79516 0 ana571_2026p002_v001 run3pp ${OUT}/${m}${SUFFIX} perfect mapscan/cm_mapscan_${m}${SUFFIX}.txt > ${OUT}/${m}/run${SUFFIX}.log 2>&1 &
  echo "started ${m} (pid $!)"
done
wait

for m in ${MAPS}; do
  echo "${m}${SUFFIX}: $(grep -c 'All done' ${OUT}/${m}/run${SUFFIX}.log) finished, $(ls ${OUT}/${m}${SUFFIX}/CMDistortionMaps_*.root 2>/dev/null | wc -l) maps"
done

# one list with everything found so far (base scan and any extensions), in config-file order
: > ${OUT}/all.list
for m in ${MAPS}; do
  for cfg in mapscan/cm_mapscan_${m}.txt mapscan/cm_mapscan_${m}_*.txt; do
    [ -f "$cfg" ] || continue
    sfx=$(basename $cfg .txt); sfx=${sfx#cm_mapscan_${m}}
    awk '!/^#/ && NF==8 {print $1}' $cfg | while read tag; do
      f=${OUT}/${m}${sfx}/CMDistortionMaps_100ADC_${tag}_g4perfect.root
      [ -f "$f" ] && echo "$f ${tag}" >> ${OUT}/all.list
    done
  done
done
root -l -b -q "PlotCMConfigProjections.C+(\"${OUT}/all.list\",\"${OUT}/mapscan_all.pdf\")" 2>&1 | grep -v "^Info" | tee ${OUT}/scores.txt
root -l -b -q "FitCMKeff.C+(\"${OUT}/all.list\",\"${OUT}/keff_fit.pdf\")" 2>&1 | grep -v "^Info" | tee ${OUT}/keff_fit.txt
