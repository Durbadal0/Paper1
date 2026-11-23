# config.sh -- settings for the simulation study (Section 4). Edit, then run the
# submit script for your scheduler. Every job reads this file.

# ---- what to run
SCENARIOS="A B C D"     # the four scenarios of Section 4
REPS=100                # simulated datasets per scenario (replicates 1..REPS)

# ---- MCMC, as stated in Section 4: 3,500 burn-in iterations followed by 6,500
#      post-burn-in iterations, posterior survival predictions kept every 10th
NUM_ITER=10000
BURN_IN=3500
THIN=10
STEP3=1                 # 1 = the 3-step algorithm (Step 3 importance weights)
NUM_TREE=200            # trees in the SBART forest (code default)

# ---- other settings: the defaults of the GitHub code
NQ=50                   # quadrature points for the Step-3 integral
FINE_N=6001             # grid for the true cumulative hazard
NB_PH=500               # bootstrap refits for the PH-frailty 95% band
NB_RSF=500              # bootstrap refits for the RSF 95% band
NT_B=250                # trees per RSF bootstrap refit
RSF_TREES=1000          # trees in the main RSF fit

# ---- resources per job (one job = one replicate of one scenario)
NCORES=8                # cores per job; the bootstrap refits are spread over them
MEM_MB=16000            # memory per job, MB (LSF: in the unit of your -M limits, usually MB)
WALLTIME="12:00"        # hours:minutes
QUEUE=""                # LSF queue or SLURM partition; empty = scheduler default
ACCOUNT=""              # SLURM account / LSF project, if your HPC requires one

# ---- R on the HPC: a command that puts Rscript on the PATH, e.g. "module load R/4.4.0".
#      If jobs report "module: command not found", start it with the module setup, e.g.
#      "source /etc/profile.d/modules.sh && module load R/4.4.0"
R_SETUP=""
# where R packages are installed if they are missing (setup.sh); empty = R's default
R_LIBS_USER_DIR="$HOME/R/sbart_sim_library"

# ---- output folders (relative to this folder)
RESULTS_DIR=results     # results/scenA/rep01.rds ...
LOG_DIR=logs

# ---- scheduler details
MAX_RUNNING=""          # at most this many jobs of an array at once (empty = no limit)
LSF_MEM_PER_CORE=0      # 1 if your LSF reads rusage[mem] per core rather than per job
