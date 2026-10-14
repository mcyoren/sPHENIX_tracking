#!/bin/bash
# Build all CM QA pdfs from the outputs under one base directory.
# usage: ./make_qa.sh [base=/sphenix/tg/tg01/hf/mitrankov/cm_configs] [outdir=<base>/qa]
MACRODIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE=${1:-/sphenix/tg/tg01/hf/mitrankov/cm_configs}
OUT=${2:-${BASE}/qa}
mkdir -p ${OUT}/lists
M=CMDistortionMaps_100ADC_

# list <name> then pairs "relative/file|label"; first entry is the reference
mklist() {
  local name=$1; shift
  : > ${OUT}/lists/${name}.list
  for e in "$@"; do
    f=${BASE}/${e%%|*}
    [ -f "$f" ] && echo "$f ${e#*|}" >> ${OUT}/lists/${name}.list || echo "  (${name}: missing ${e%%|*})"
  done
}
qa() {  # qa <name> <title>
  [ -s ${OUT}/lists/$1.list ] || return
  echo "== $1"
  root -l -b -q "${MACRODIR}/CMQA.C+(\"${OUT}/lists/$1.list\",\"${OUT}/$1.pdf\",\"$2\",0)" > ${OUT}/lists/$1.txt 2>&1
  grep -E "S: <dR>" ${OUT}/lists/$1.txt | cut -c1-200
}

mklist QA_01_survey_noGarfield \
  "g4perfect/${M}noGarfield_g4perfect.root|no Garfield, G4 perfect" \
  "g4survey/${M}noGarfield_g4survey.root|no Garfield, G4 survey, global" \
  "g4survey_local/${M}noGarfield_g4survey_local.root|no Garfield, G4 survey, local"
qa QA_01_survey_noGarfield "Survey geometry without Garfield: global vs local coordinates"

mklist QA_02_survey_Garfield \
  "g4perfect/${M}cont_cdb_g4perfect.root|G4 perfect" \
  "g4survey/${M}cont_cdb_g4survey.root|G4 survey, global, Garf container" \
  "g4survey_local/${M}cont_cdb_g4survey_local.root|G4 survey, local, Garf container" \
  "g4survey/${M}hard_cdb_g4survey.root|G4 survey, global, Garf hard" \
  "g4survey_local/${M}hard_cdb_g4survey_local.root|G4 survey, local, Garf hard" \
  "g4survey/${M}perf_cdb_g4survey.root|G4 survey, global, Garf perfect" \
  "g4survey_local/${M}perf_cdb_g4survey_local.root|G4 survey, local, Garf perfect"
qa QA_02_survey_Garfield "Survey geometry with Garfield (CDB 2D map, CDB keff)"

mklist QA_03_Garfield_tpcgeom \
  "g4perfect/${M}cont_cdb_g4perfect.root|G4 perfect, Garf container" \
  "g4perfect/${M}hard_cdb_g4perfect.root|G4 perfect, Garf hard" \
  "g4perfect/${M}perf_cdb_g4perfect.root|G4 perfect, Garf perfect"
qa QA_03_Garfield_tpcgeom "PHGarfield TPC placement mode, G4 perfect"

T=tiltscan/${M}cont_cdb_g4tilt_-0p298_
mklist QA_04_tilt_ry "${T}0_local.root|rot_y 0" "${T}0p37_local.root|rot_y 0.37" "${T}0p74_local.root|rot_y 0.74" \
  "${T}1p11_local.root|rot_y 1.11" "${T}1p48_local.root|rot_y 1.48 (survey)"
qa QA_04_tilt_ry "TPC tilt scan rot_y [mrad], rot_x = -0.298, local coordinates, Garfield container"

T=tiltscan_rx/${M}cont_cdb_g4tilt_
mklist QA_05_tilt_rx "${T}-0p298_0p85_local.root|rot_x -0.298 (survey)" "${T}-0p6_0p85_local.root|rot_x -0.6" \
  "${T}0_0p85_local.root|rot_x 0" "${T}0p16_0p85_local.root|rot_x 0.16" "${T}0p3_0p85_local.root|rot_x 0.3" "${T}0p6_0p85_local.root|rot_x 0.6"
qa QA_05_tilt_rx "TPC tilt scan rot_x [mrad], rot_y = 0.85, local coordinates, Garfield container"

T=keffscan_perfect/${M}
mklist QA_06_keff_absolute "${T}keff_cdb_g4perfect.root|keff CDB (1.39/1.38)" "${T}keff0p0_g4perfect.root|keff 0" \
  "${T}keff0p5_g4perfect.root|keff 0.5" "${T}keff1p0_g4perfect.root|keff 1.0" "${T}keff1p5_g4perfect.root|keff 1.5" \
  "${T}keff2p0_g4perfect.root|keff 2.0" "${T}keff3p0_g4perfect.root|keff 3.0"
qa QA_06_keff_absolute "Absolute keff scan, CDB 2D map, G4 perfect"

KPTS="0p900 0p950 0p975 1p000 1p025 1p050 1p100 1p150 1p200 1p300 1p400"
MAPS="cdb2d alpha1p25 alpha1p6 alpha2p0"

# relative path of one scan point: <dir>/<map>[_ext]/<M><map>_k<k>_<suffix>.root
scanfile() {  # scanfile <dir> <suffix> <map> <k>
  for sub in $3 $3_ext; do
    f=$1/${sub}/${M}$3_k$4_$2.root
    [ -f ${BASE}/$f ] && { echo $f; return; }
  done
  echo $1/$3/${M}$3_k$4_$2.root
}

# all QA for one map x keff scan set
scanset() {  # scanset <set name> <dir> <file suffix> <description>
  local set=$1 dir=$2 sfx=$3 desc=$4
  [ -d ${BASE}/${dir} ] || { echo "  (no ${dir}, skipping ${set})"; return; }
  for m in ${MAPS}; do
    local args=("$(scanfile ${dir} ${sfx} ${m} 1p000)|${m} k x1.000")
    for k in ${KPTS}; do
      [ $k = 1p000 ] || args+=("$(scanfile ${dir} ${sfx} ${m} ${k})|${m} k x${k/p/.}")
    done
    mklist QA_07_keff_${m}_${set} "${args[@]}"
    qa QA_07_keff_${m}_${set} "Relative keff scan (x CDB keff), map ${m}, ${desc}"
  done

  local args=()
  for m in ${MAPS}; do args+=("$(scanfile ${dir} ${sfx} ${m} 1p000)|${m}"); done
  args+=("${dir}/cdb2d/${M}noGarfield_${sfx}.root|noGarfield")
  mklist QA_08_maps_at_kCDB_${set} "${args[@]}"
  qa QA_08_maps_at_kCDB_${set} "Space-charge maps at the CDB keff, ${desc}"

  : > ${OUT}/lists/keff_fit_${set}.list
  for m in ${MAPS}; do
    for k in ${KPTS}; do
      f=${BASE}/$(scanfile ${dir} ${sfx} ${m} ${k})
      [ -f "$f" ] && echo "$f ${m}_k${k}" >> ${OUT}/lists/keff_fit_${set}.list
    done
  done
  echo "== QA_09_keff_fit_${set}"
  root -l -b -q "${MACRODIR}/FitCMKeff.C+(\"${OUT}/lists/keff_fit_${set}.list\",\"${OUT}/QA_09_keff_fit_${set}.pdf\")" 2>&1 | grep -v "^Info\|^$\|Processing" | tee ${OUT}/QA_09_keff_fit_${set}.txt

  # maps at their fitted keff (nearest scan point, M1 fit averaged over sides)
  args=()
  for m in ${MAPS}; do
    kbest=$(awk -v m=$m '$1==m && ($2=="South" || $2=="North") {s+=$7; n++} END {if (n) printf "%.4f", s/n}' ${OUT}/QA_09_keff_fit_${set}.txt)
    [ -n "$kbest" ] || continue
    best=$(for k in ${KPTS}; do echo "$k $(awk -v a=${k/p/.} -v b=$kbest 'BEGIN{d=a-b; print (d<0?-d:d)}')"; done | sort -k2 -g | head -1 | cut -d' ' -f1)
    args+=("$(scanfile ${dir} ${sfx} ${m} ${best})|${m} k x${best/p/.} (fit ${kbest})")
  done
  mklist QA_10_maps_at_fitted_keff_${set} "${args[@]}"
  qa QA_10_maps_at_fitted_keff_${set} "Space-charge maps at their fitted keff (nearest scan point), ${desc}"
}

scanset ideal   mapscan                       g4perfect                 "ideal TPC, no tilt w.r.t. B"
scanset tiltfit mapscan_tilt_0p16_0p85_local  g4tilt_0p16_0p85_local    "fitted TPC tilt (rx 0.16, ry 0.85 mrad), local coordinates"
scanset survey  mapscan_survey_local          g4survey_local            "survey TPC tilt (rx -0.30, ry 1.48 mrad), local coordinates"

# same map and keff, three TPC geometries
for m in cdb2d alpha1p6; do
  mklist QA_12_geometry_${m}_kCDB \
    "$(scanfile mapscan g4perfect ${m} 1p000)|${m} k x1.000, ideal TPC" \
    "$(scanfile mapscan_tilt_0p16_0p85_local g4tilt_0p16_0p85_local ${m} 1p000)|${m} k x1.000, fitted tilt" \
    "$(scanfile mapscan_survey_local g4survey_local ${m} 1p000)|${m} k x1.000, survey tilt"
  qa QA_12_geometry_${m}_kCDB "TPC geometry seen by Garfield: ideal vs fitted vs survey tilt, map ${m}, CDB keff"
done

mklist QA_11_tilt_vs_ideal \
  "g4perfect/${M}cont_cdb_g4perfect.root|ideal TPC (no tilt w.r.t. B)" \
  "tiltscan_rx/${M}cont_cdb_g4tilt_0p16_0p85_local.root|fitted tilt rx 0.16 ry 0.85 mrad, local" \
  "g4survey_local/${M}cont_cdb_g4survey_local.root|survey tilt rx -0.30 ry 1.48 mrad, local"
qa QA_11_tilt_vs_ideal "Garfield with ideal TPC vs tilted TPC (CDB 2D map, CDB keff, local coordinates)"
echo "pdfs in ${OUT}"
