#!/bin/bash
# collect_results.sh -- summarise the finished replicates: one summary per scenario
# (results/scenX/RESULTS.txt and figures) and the outputs for the paper
# (paper_outputs/). Can be run at any time; it uses the replicates finished so far.
set -eo pipefail
cd "$(dirname "$0")"
source ./config.sh
[ -n "$R_SETUP" ] && eval "$R_SETUP"
[ -n "$R_LIBS_USER_DIR" ] && export R_LIBS="$R_LIBS_USER_DIR${R_LIBS:+:$R_LIBS}"
for S in $SCENARIOS; do
  d="$RESULTS_DIR/scen$S"
  ls "$d"/rep*.rds > /dev/null 2>&1 || { echo "scenario $S: no results yet"; continue; }
  if Rscript code/sim_aggregate.R "$d" > "$d/RESULTS.txt" 2>&1; then
    echo "scenario $S: $d/RESULTS.txt"
  else
    echo "scenario $S: sim_aggregate.R FAILED; the end of $d/RESULTS.txt:"; tail -n 5 "$d/RESULTS.txt"
  fi
done
Rscript code/make_paper_outputs.R "$RESULTS_DIR" paper_outputs
