#!/bin/bash
# job.sh SCENARIO [REP] -- run one replicate of one scenario.
# Called by submit_lsf.sh, submit_slurm.sh and run_one_machine.sh. REP defaults to
# the array index of the scheduler (LSB_JOBINDEX or SLURM_ARRAY_TASK_ID).
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
source ./config.sh
SCEN=$1
REP=${2:-${LSB_JOBINDEX:-${SLURM_ARRAY_TASK_ID:-}}}
[ -n "$SCEN" ] && [ -n "$REP" ] || { echo "usage: job.sh SCENARIO REP"; exit 1; }
[ -n "$R_SETUP" ] && eval "$R_SETUP"
[ -n "$R_LIBS_USER_DIR" ] && export R_LIBS="$R_LIBS_USER_DIR${R_LIBS:+:$R_LIBS}"
# one thread for BLAS; NCORES for ranger and for the bootstrap refits
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export R_RANGER_NUM_THREADS=$NCORES NCORES
OUT="$RESULTS_DIR/scen$SCEN"
mkdir -p "$OUT"
RDS=$(printf "%s/rep%02d.rds" "$OUT" "$REP")
if [ -f "$RDS" ]; then echo "$RDS exists: nothing to do"; exit 0; fi
echo "scenario $SCEN, replicate $REP, $NCORES cores, host $(hostname), start $(date)"
Rscript code/sim_worker.R "$REP" "$OUT" "$NUM_ITER" "$BURN_IN" "$NUM_TREE" "$SCEN" "$STEP3" \
  "$THIN" "$NQ" "$FINE_N" "$NB_PH" "$NB_RSF" "$NT_B" "$RSF_TREES" "$NCORES"
echo "end $(date)"
