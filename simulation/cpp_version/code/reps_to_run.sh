#!/bin/bash
# reps_to_run.sh SCENARIO [--missing] -- the replicate numbers to submit, as an
# array index list (1-100, or 3,17,42 for the missing ones only)
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source ./config.sh
S=$1
if [ "$2" != "--missing" ]; then echo "1-$REPS"; exit 0; fi
miss=""
for r in $(seq 1 "$REPS"); do
  [ -f "$(printf "%s/scen%s/rep%02d.rds" "$RESULTS_DIR" "$S" "$r")" ] || miss="$miss,$r"
done
echo "${miss#,}"
