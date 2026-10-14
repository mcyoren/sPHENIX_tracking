#!/bin/bash
# Ben's BenLaserAggregatedClustering_simple.sh with own install and config arguments
# usage: CMLaserAggregatedClustering_configs.sh run segStart prod type outdir [g4geom=cdb|survey|perfect] [configfile] [doQA=0|1] [useGlobal=0|1, default 0 = local]

export USER="$(id -u -n)"
export LOGNAME=${USER}
export MYINSTALL=${MYINSTALL:-/sphenix/tg/tg01/hf/mitrankov/code/install_cm}
MACRODIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "Done with setup"

runnumber=$1
segmentStart=$2
prod=$3
type=$4
outdir=$5
g4geom=${6:-cdb}
configfile=${7:-${MACRODIR}/cm_garfield_configs.txt}
doQA=${8:-0}
if [ "$doQA" = "1" ]; then doQA=true; else doQA=false; fi
useGlobal=${9:-0}
if [ "$useGlobal" = "1" ]; then useGlobal=true; else useGlobal=false; fi

mkdir -p ${outdir}
cd ${MACRODIR}

outfile="DST_LASER_AGGREGATED_CLUSTER_${type}_${prod}-${runnumber}-${segmentStart}_g4${g4geom}.root"

filelist=${outdir}/run${runnumber}_seg${segmentStart}_g4${g4geom}_clusterList.list
rm -f $filelist

for i in `seq ${segmentStart} $((segmentStart+99))`; do
    ls /sphenix/tg/tg01/jets/bkimelman/BenProduction/Oct9_2026_LaserAggregatedPads/DST_LASER_AGGREGATED_PADS_${type}_${prod}-000${runnumber}-$(printf "%05d" $i).root >> $filelist
done

root -b -q Fun4All_LaserAggregatedClustering_configs.C\(\"${filelist}\",\"${outdir}\",\"${outfile}\",\"${g4geom}\",\"${configfile}\",${doQA},${useGlobal}\)

rm $filelist
