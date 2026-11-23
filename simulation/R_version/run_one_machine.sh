#!/bin/bash
# run_one_machine.sh [JOBS_AT_ONCE] -- run all replicates on one machine without a
# scheduler (a large server, or an interactive HPC node). Each replicate uses NCORES
# cores; JOBS_AT_ONCE replicates run side by side (default: cores / NCORES).
# Replicates that already have a results file are skipped, so it can be restarted.
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
NC=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)
J=${1:-$(( NC / NCORES > 0 ? NC / NCORES : 1 ))}
echo "running with $J replicates at once, $NCORES cores each"
for S in $SCENARIOS; do for r in $(seq 1 "$REPS"); do echo "$S $r"; done; done |
  xargs -P "$J" -L 1 sh -c 'bash code/job.sh "$0" "$1" > "'"$LOG_DIR"'/${0}_${1}.out" 2>&1'
echo "all done; now run: bash collect_results.sh"
