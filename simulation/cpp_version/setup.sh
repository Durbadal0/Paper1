#!/bin/bash
# setup.sh -- run once on the HPC (login node or a short interactive job):
#   1. installs missing R packages (into $R_LIBS_USER_DIR);
#   2. compiles the C++ engine once, into code/cpp/build (the jobs only load it);
#   3. checks the engine: truncated-normal draws, the forest code, and one short
#      replicate run by both the R and the C++ version, compared field by field.
set -eo pipefail
cd "$(dirname "$0")"
source ./config.sh
[ -n "$R_SETUP" ] && eval "$R_SETUP"
if [ -n "$R_LIBS_USER_DIR" ]; then mkdir -p "$R_LIBS_USER_DIR"; export R_LIBS="$R_LIBS_USER_DIR${R_LIBS:+:$R_LIBS}"; fi
Rscript code/install_packages.R
echo "compiling the C++ engine ..."
Rscript -e 'source("code/cpp/load_engine.R"); load_engine("code", rebuild = TRUE); cat("engine compiled\n")'
Rscript verify/check_truncnorm.R
Rscript verify/check_forest_identity.R 2>&1 | grep -v "Finishing iteration"
NCORES=2 bash verify/run_comparison.sh quick
echo "setup OK. Submit with: bash submit_lsf.sh   (or submit_slurm.sh, or run_one_machine.sh)"
