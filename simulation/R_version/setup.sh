#!/bin/bash
# setup.sh -- run once on the HPC (login node is fine): installs missing R packages
# and runs one very short replicate as a test (about a minute).
set -eo pipefail
cd "$(dirname "$0")"
source ./config.sh
[ -n "$R_SETUP" ] && eval "$R_SETUP"
if [ -n "$R_LIBS_USER_DIR" ]; then mkdir -p "$R_LIBS_USER_DIR"; export R_LIBS="$R_LIBS_USER_DIR${R_LIBS:+:$R_LIBS}"; fi
Rscript code/install_packages.R
T=$(mktemp -d)
echo "test replicate (short settings) in $T ..."
OMP_NUM_THREADS=1 R_RANGER_NUM_THREADS=2 Rscript code/sim_worker.R 1 "$T" 300 200 20 A 1 10 10 1001 6 6 30 60 2 > "$T/test.log" 2>&1 \
  || { echo "TEST FAILED, see $T/test.log"; exit 1; }
grep "DONE" "$T"/rep01.log
echo "setup OK. Submit with: bash submit_lsf.sh   (or submit_slurm.sh, or run_one_machine.sh)"
