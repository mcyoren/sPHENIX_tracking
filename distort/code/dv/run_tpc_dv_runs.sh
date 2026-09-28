#!/bin/bash
set -u

# Run PHGarfield drift-velocity scans for a predefined set of runs.
#
# Usage:
#   ./run_tpc_dv_runs.sh short
#   ./run_tpc_dv_runs.sh long
#
# Optional environment variables:
#   MACRO=/full/path/to/Fun4All_TpcDriftVelocityScan.C
#   ROOTDIR=/sphenix/tg/tg01/hf/mitrankov/dv/rootfiles
#   LOGDIR=/sphenix/tg/tg01/hf/mitrankov/dv/logs
#   FORCE=1     # rerun even if output ROOT file already exists
#
# The scan macro is called with only the run number.  This assumes your
# macro's default event input is sufficient to drive se->run(1), while the
# run number / Tracking_Geometry / PHGarfield CDB inputs are selected by CDB.

MODE="${1:-short}"

MACRO="${MACRO:-Fun4All_TpcDriftVelocityScan.C}"
ROOTDIR="${ROOTDIR:-/sphenix/tg/tg01/hf/mitrankov/dv/rootfiles}"
LOGDIR="${LOGDIR:-/sphenix/tg/tg01/hf/mitrankov/dv/logs}"
FORCE="${FORCE:-0}"

mkdir -p "${ROOTDIR}" "${LOGDIR}"

SHORT_RUNS=(
  79507 79508 79509 79510 79511 79512 79513 79514 79515 79516
  79522 79523 79524 79525 79526 79527 79528 79529 79530
)

LONG_RUNS=(
  67777 67778 67779 67780 67781 67782 67784 67785 67786
  72309 72310 72311 72312 72313 72314 72316
  74502 74503 74504 74505 74506 74507 74509 74510
  75570 75574 76905 75405 75391
  78014 78015 78016 78017 78018 78019 78020 78029 78030 78032 78033 78034 78035 78036
  79506 79507 79508 79509 79510 79511 79512 79513 79514 79515 79516
  79522 79523 79524 79525 79526 79527 79528 79529 79530
  79563 79564 79565 79566 79567 79568 79569 79570 79571 79572
  80858 80859 80860 80861 80863 80864 80865 80866 80867
  80878 80879 80880 80881
  80891 80892 80893 80895
  80901 80902 80903 80904 80905 80906 80908 80909 80910 80911 80912 80913 80914 80915
  81556 81558 81559 81560 81561 81562 81563 81564 81565 81566
  81575 81576 81577 81580 81581 81583 81584 81585 81586
  81600 81601 81602 81603 81605 81607 81608 81609 81610 81611 81612
  82391 82420 82435 82436 82439 82440 82441 82442 82443
  82624 82625 82626 82627 82628 82629 82630 82631 82632
)

case "${MODE}" in
  short)
    RUNS=("${SHORT_RUNS[@]}")
    ;;
  long)
    RUNS=("${LONG_RUNS[@]}")
    ;;
  *)
    echo "Usage: $0 {short|long}"
    exit 2
    ;;
esac

echo "Mode      : ${MODE}"
echo "Macro     : ${MACRO}"
echo "ROOT dir  : ${ROOTDIR}"
echo "Log dir   : ${LOGDIR}"
echo "N runs    : ${#RUNS[@]}"
echo

n_ok=0
n_fail=0
n_skip=0

for run in "${RUNS[@]}"; do
  outfile="${ROOTDIR}/tpc_dv_scan_run${run}.root"
  logfile="${LOGDIR}/tpc_dv_scan_run${run}.log"

  if [[ "${FORCE}" != "1" && -s "${outfile}" ]]; then
    echo "[SKIP] run ${run}: ${outfile} already exists"
    ((n_skip+=1))
    continue
  fi

  echo "[RUN ] ${run}"
  echo "       log: ${logfile}"

  # Run in a separate ROOT process for each run so CDB/PHGarfield state is reset.
  root --web=off -l -b -q "${MACRO}(${run})" >"${logfile}" 2>&1
  status=$?

  if [[ ${status} -eq 0 && -s "${outfile}" ]]; then
    echo "[ OK ] ${run}"
    ((n_ok+=1))
  else
    echo "[FAIL] ${run} (ROOT exit=${status})"
    echo "       see ${logfile}"
    ((n_fail+=1))
  fi
done

echo
echo "Finished."
echo "  OK      : ${n_ok}"
echo "  skipped : ${n_skip}"
echo "  failed  : ${n_fail}"

if (( n_fail > 0 )); then
  exit 1
fi
