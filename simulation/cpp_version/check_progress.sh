#!/bin/bash
# check_progress.sh -- how many replicates are finished, and how long they took
cd "$(dirname "$0")"
source ./config.sh
for S in $SCENARIOS; do
  d="$RESULTS_DIR/scen$S"
  n=$(ls "$d"/rep*.rds 2>/dev/null | wc -l | tr -d ' ')
  mins=$(grep -h "DONE" "$d"/rep*.log 2>/dev/null | sed -E 's/.*\| ([0-9.]+) min.*/\1/' |
         awk '{s+=$1; if($1>m)m=$1} END {if (NR) printf "mean %.0f min, longest %.0f min", s/NR, m}')
  running=$(for f in "$d"/rep*.log; do [ -f "$f" ] || continue; b=${f%.log}; [ -f "$b.rds" ] || echo "$f"; done | wc -l | tr -d ' ')
  echo "scenario $S: $n of $REPS finished ($mins); started but not finished: $running"
done
echo "Failed jobs leave an error message in $LOG_DIR/*.err. Once no job of an array is"
echo "queued or running any more (bjobs / squeue), resubmit what is missing with --missing."
