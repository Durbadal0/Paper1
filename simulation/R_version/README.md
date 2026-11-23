# Simulation study (Section 4): R version

This folder runs the four simulation scenarios of Section 4, with 100 simulated datasets each. Each dataset is one job, so the 400 jobs can run in parallel. When they have finished, one command writes Figure 2, Figures S1–S4, Table S1 and the numbers of Section 4.

The folder is self-contained: copy it to the cluster, edit `config.sh`, and run the scripts from inside the folder.

## Steps

```bash
# 1. once: R packages (installed into $R_LIBS_USER_DIR if missing) and a one-minute test
bash setup.sh

# 2. submit the 400 jobs: 4 job arrays of 100 (pick your scheduler)
bash submit_lsf.sh        # LSF (bsub)
bash submit_slurm.sh      # SLURM (sbatch)
bash run_one_machine.sh   # no scheduler: one large server or an interactive node

# 3. while they run
bash check_progress.sh

# 4. if some jobs failed or ran out of time: submit only the missing replicates
bash submit_lsf.sh --missing      # or: bash submit_slurm.sh --missing

# 5. when they have finished (or at any time, with the replicates done so far)
bash collect_results.sh
```

Before step 1, edit `config.sh`:

- `R_SETUP`: how R is loaded on your cluster, e.g. `module load R/4.4.0`.
- `QUEUE`, and `ACCOUNT` if your cluster needs one.
- `NCORES`, `MEM_MB`, `WALLTIME` if the defaults do not suit your cluster.

## What one job does

`code/sim_worker.R` runs one replicate of one scenario:

- **Data.** It simulates one dataset with seed 20260814 + replicate number.
- **Competing methods.** It fits the RSF and the PH-frailty model.
- **SBART.** It runs the three steps of the method:
  - Step 1: CAR model for M;
  - Step 2: SBART survival model;
  - Step 3: importance weights.
- **Accuracy.** It computes the AMSE and the AES curves.
- **Intervals.** It computes each method's 95% pointwise interval: the SBART credible band, and bootstrap bands from 500 refits for PH-frailty and RSF.
- **Output.** It saves everything in `results/scenX/repNN.rds`, with a log in `results/scenX/repNN.log`.

The script is the single-core script `scenA.R` (one file now serves all four scenarios). It differs only in how it runs, not in what it computes:

- **Bootstrap on several cores.** The 1,000 bootstrap refits (500 PH-frailty, 500 RSF) run on `NCORES` cores instead of one. The resamples are still drawn in the original order. We checked that the output file is identical to that of a one-core run of `scenA.R` and `scenC.R`.
- **Safe saving.** The result file is first written as `repNN.rds.part` and renamed when complete. A job killed while saving therefore leaves no half-written file that would count as finished.

## Settings (`config.sh`)

| Setting | Value |
|---|---|
| Replicates per scenario | 100 |
| Iterations, burn-in | 10,000, 3,500 (6,500 kept) |
| Survival predictions kept | every 10th iteration (650 draws) |
| Algorithm | 3 steps (`STEP3=1`) |
| Trees K | 200 |
| Bootstrap refits PH / RSF | 500 / 500 (250 trees each) |
| RSF trees, main fit | 1,000 |
| Step-3 quadrature points, fine grid | 50, 6,001 |

## Time and memory per job

These figures were measured on a Mac (Apple silicon, R 4.6.1) for one replicate at the full settings (10,000 iterations, K = 200, 650 kept draws), with four such jobs running at once.

| | Scenario A | Scenario D |
|---|---|---|
| Steps 1–3 of the SBART method (one core; this part cannot be split) | 44 min | 41 min |
| Peak memory | 8.2 GB | 5.1 GB |

In addition:

- **Bootstrap refits.** The 500 PH-frailty refits take about 2 min of core time. The 500 RSF refits take about 2 hours of core time, divided over `NCORES`.
- **Cluster cores.** Cluster cores are often slower than this. `WALLTIME=12:00` and `MEM_MB=16000` leave a margin. With the C++ version, the median time per replicate on the cluster used for the paper was 230–353 minutes (`../outputs/section4_numbers.txt`).

Each job reserves its `NCORES` cores for the whole run, and the SBART part uses only one of them.

If your cluster limits how many jobs one user can run, set `MAX_RUNNING`. The limit applies to each of the four arrays, so up to 4 × `MAX_RUNNING` jobs run at once.

## Outputs (`collect_results.sh`)

For each scenario, `results/scenX/` gets:

- `RESULTS.txt`: AMSE summary, intervals and importance-sampling ESS;
- `amse_by_rep.csv`;
- `AES_X.pdf` and `AMSE_X.pdf`;
- `intervals_by_rep.rds`.

These come from `code/sim_aggregate.R`, the per-scenario aggregator `scenA_aggregate.R`. It skips the PNG copies of the figures when R cannot write PNG files; the PDFs are always written.

For the paper, `paper_outputs/` gets the following (the copies used in the paper are in `../outputs/`):

| File | In the paper |
|---|---|
| `Figure2_AMSE_boxplots.pdf` | Figure 2: per-dataset AMSE, one panel per scenario |
| `FigureS1_AES_A.pdf` … `FigureS4_AES_D.pdf` | Figures S1–S4; band = 2.5%–97.5% range of the replicate AES curves |
| `FigureS*_AES_*_ownband.pdf` | not in the paper; band = each method's own 95% interval, averaged over replicates |
| `TableS1_coverage.csv`, `.tex` | Table S1: coverage of each method's own 95% pointwise interval |
| `amse_summary.csv`, `amse_by_replicate.csv` | AMSE mean, SD, median and quartiles; the AMSE of every replicate |
| `section4_numbers.txt` | the numbers of Section 4, with replicate counts, censoring, settings and run times |

## Folder contents

```
README.md                 this file
config.sh                 all settings
setup.sh                  package check/installation and a short test run
submit_lsf.sh             LSF job arrays (bsub)
submit_slurm.sh           SLURM job arrays (sbatch)
run_one_machine.sh        all jobs on one machine, without a scheduler
check_progress.sh         finished replicates and run times
collect_results.sh        summaries and the paper outputs
code/sim_worker.R         one replicate
code/sim_aggregate.R      per-scenario summary
code/make_paper_outputs.R Figure 2, Figures S1-S4, Table S1, Section 4 numbers
code/job.sh               one job (called by the submit scripts)
code/reps_to_run.sh       which replicates to submit
code/install_packages.R   installs SoftBart, MASS, truncnorm, survival, ranger if missing
```

The checks on the Mac used R 4.6.1, SoftBart 1.0.3, truncnorm 1.0.9, ranger 0.18.0, survival 3.8.6 and MASS 7.3.65. With other versions or on other machines the code runs the same way, but the later digits of the results can differ (see `../README.md`).
