#!/bin/bash
# run_comparison.sh [quick|small|full] -- runs the R version (verify/sim_worker_R.R,
# the same file as sim_worker.R of the R folder) and the C++ version on the same
# replicates and compares every saved quantity.
#   quick : scenario A, replicate 1, short chains (about 2 minutes)
#   small : scenarios A-D, replicates 1-3, short chains (about 15 minutes on 4 cores)
#   full  : scenario A, replicate 1, the settings of config.sh but 40/16 bootstrap
#           refits (about 1 hour)
set -eo pipefail
cd "$(dirname "$0")/.."
NC_ARG=${NCORES:-}          # NCORES from the environment, else 4 (config.sh's NCORES is not used here)
source ./config.sh
[ -n "$R_SETUP" ] && eval "$R_SETUP"
[ -n "$R_LIBS_USER_DIR" ] && export R_LIBS="$R_LIBS_USER_DIR${R_LIBS:+:$R_LIBS}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
MODE=${1:-quick}; NC=${NC_ARG:-4}; export R_RANGER_NUM_THREADS=$NC
case $MODE in
  quick) SC="A"; RP="1";     SET="300 200 20 SC 1 10 10 1001 6 6 30 60" ;;
  small) SC="A B C D"; RP="1 2 3"; SET="600 300 50 SC 1 10 20 2001 10 6 50 100" ;;
  full)  SC="A"; RP="1";     SET="$NUM_ITER $BURN_IN $NUM_TREE SC $STEP3 $THIN $NQ $FINE_N 40 16 $NT_B $RSF_TREES" ;;
  *) echo "usage: run_comparison.sh [quick|small|full]"; exit 1 ;;
esac
OUT=verify/out_$MODE; rm -rf "$OUT"; mkdir -p "$OUT"
for S in $SC; do for r in $RP; do
  A=${SET/SC/$S}
  echo "scenario $S replicate $r: R version ..."
  Rscript verify/sim_worker_R.R "$r" "$OUT/R/scen$S" $A "$NC" > "$OUT/R_${S}_$r.out" 2>&1
  echo "scenario $S replicate $r: C++ version ..."
  Rscript code/sim_worker_cpp.R "$r" "$OUT/cpp/scen$S" $A "$NC" > "$OUT/cpp_${S}_$r.out" 2>&1
done; done
for S in $SC; do
  echo; echo "== scenario $S"
  Rscript verify/compare_outputs.R "$OUT/R/scen$S" "$OUT/cpp/scen$S" | tee "$OUT/compare_$S.txt"
done
