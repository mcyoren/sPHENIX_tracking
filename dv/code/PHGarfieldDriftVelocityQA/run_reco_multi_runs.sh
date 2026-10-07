#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}



#source /opt/sphenix/core/bin/sphenix_setup.sh new
#source /sphenix/u/mitrankova/.login_TPC_SA
source /sphenix/u/mitrankov/.login_PR_3D
#printenv

nEvents="$1"
#"ana537_nocdbtag_v001"

# Run number is supplied by the Condor submit file.
runnumber="$2"

# Set collision system and CDB tag based on run number
case "${runnumber}" in

  79507|79508|79509|79510|79511|79512|79513|79514|79515|79516|79526|79528|79529|79523|81558|80892|80909|81566)
    collision="run3pp"
    cdbtag="ana532_nocdbtag_v001"
    ;;

  75391|75405|76905)
    collision="run3auau"
    cdbtag="ana514_nocdbtag_v001"
    ;;

  82525)
    collision="run3oo"
    cdbtag="ana537_nocdbtag_v001"
    ;;

  *)
    echo "ERROR: Unknown run number ${runnumber}"
    echo "No collision system / CDB tag configured for this run."
    exit 1
    ;;
esac

#collision="run3pp"
#cdbtag="ana532_nocdbtag_v001"
OutDir=/sphenix/tg/tg01/hf/mitrankov/ibf/output/
name="DriftVelocityQA"

#runnumber=76905
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905_new_8_380_no_mv_rot_Btest
#name="HITS_6x6_1mrad_AuAu"

#runnumber=82626
#collision="run3oo"
#cdbtag="ana537_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/82626_00
#name="HITS_OO"

#runnumber=75405
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/75405_new
#name="HITS_AuAuFieldOn"

#runnumber=75391
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/75391_new
#name="HITS_AuAuFieldOn"

segment="$3"
nSkip="$4"

echo ========================================================
echo Running Fun4All_raw_hit_TPC_DriftVelocityQA.C
echo "  nEvents             = ${nEvents}"
echo "  nSkip               = ${nSkip}"
echo "  Run Number          = ${runnumber}"
echo "  Segment             = ${segment}"
echo "  OutDir              = ${OutDir}"
echo "  collision           = ${collision}"
echo "  cdbtag              = ${cdbtag}"
echo "  name                = ${name}"
echo ========================================================


 mkdir -p "${OutDir}/QA"
 root.exe -l -b << EOF
        .x Fun4All_raw_hit_TPC_DriftVelocityQA.C(${nEvents},${runnumber}, ${segment},"${OutDir}",${nSkip},"${collision}","${cdbtag}","${name}")
EOF

status=$?
exit "$status"
