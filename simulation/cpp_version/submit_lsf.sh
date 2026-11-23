#!/bin/bash
# submit_lsf.sh [--missing] -- submit one LSF job array per scenario (bsub).
# --missing submits only the replicates that have no results file yet.
set -eo pipefail
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
MEMREQ=$MEM_MB
[ "$LSF_MEM_PER_CORE" = "1" ] && MEMREQ=$(( (MEM_MB + NCORES - 1) / NCORES ))
for S in $SCENARIOS; do
  IDX=$(bash code/reps_to_run.sh "$S" "$1")
  if [ -z "$IDX" ]; then echo "scenario $S: all $REPS replicates done"; continue; fi
  bsub -J "sim${S}[${IDX}]${MAX_RUNNING:+%$MAX_RUNNING}" \
       -n "$NCORES" -R "span[hosts=1]" -R "rusage[mem=$MEMREQ]" -M "$MEM_MB" \
       -W "$WALLTIME" ${QUEUE:+-q "$QUEUE"} ${ACCOUNT:+-P "$ACCOUNT"} \
       -o "$LOG_DIR/${S}_%I.out" -e "$LOG_DIR/${S}_%I.err" \
       bash code/job.sh "$S"
done
echo "Check with: bjobs    and    bash check_progress.sh"
