#!/bin/bash
# submit_slurm.sh [--missing] -- submit one SLURM job array per scenario (sbatch).
# --missing submits only the replicates that have no results file yet.
set -eo pipefail
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
for S in $SCENARIOS; do
  IDX=$(bash code/reps_to_run.sh "$S" "$1")
  if [ -z "$IDX" ]; then echo "scenario $S: all $REPS replicates done"; continue; fi
  sbatch --job-name="sim$S" --array="${IDX}${MAX_RUNNING:+%$MAX_RUNNING}" \
         --nodes=1 --ntasks=1 --cpus-per-task="$NCORES" --mem="${MEM_MB}M" \
         --time="${WALLTIME}:00" ${QUEUE:+--partition="$QUEUE"} ${ACCOUNT:+--account="$ACCOUNT"} \
         --output="$LOG_DIR/${S}_%a.out" --error="$LOG_DIR/${S}_%a.err" \
         --wrap="bash code/job.sh $S"
done
echo "Check with: squeue -u \$USER    and    bash check_progress.sh"
