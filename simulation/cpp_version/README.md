# Simulation study (Section 4): C++ version

This folder runs the same simulation as the R version (`../R_version`): four scenarios × 100 datasets, one job per dataset. The SBART sampler itself runs in C++. The outputs of the paper (`../outputs/`) were produced with this version.

**The C++ version gives the same results as the R version.** It uses the same seeds, the same data and the same random draws. `VERIFICATION.md` reports the checks, run on a Mac, in which the two versions gave identical results. `setup.sh` repeats the main checks on the machine where it is run.

## What is in C++ and what stays in R

| Part | Where |
|---|---|
| Data generation, true survival curves | R (same code as the R version) |
| RSF (`ranger`) and PH-frailty (`survival::coxph`) fits and their bootstrap bands | R (these packages are compiled code already) |
| Step 1: CAR model for M (12,000 iterations) | **C++**, `sbart_step1()` |
| Steps 2–3: SBART sampler, Poisson-thinning augmentation, λ0, σ1², ρ1, frailties, importance weights, AMSE and AES running sums | **C++**, `sbart_step2()` |
| AMSE, AES and interval summaries, saving | R (same code as the R version) |

`code/sim_worker_cpp.R` is `sim_worker.R` of the R version with only two loops replaced by calls to the C++ engine (`code/cpp/sbart_engine.cpp`). These are the Step-1 loop and the Step-2/Step-3 loop.

The engine uses the SoftBart forest code itself (`code/cpp/softbart/`). It draws from R's random-number generator, in the same order as the R loops. The large prediction blocks (the 15 × 512 × 101 AMSE grid, the AES grid and the Step-3 quadrature) are spread over `NCORES` threads. Every row is computed exactly as SoftBart computes it, so the number of threads does not change any result.

## Steps

```bash
# 1. once: R packages, compile the C++ engine, check it against the R version (~5 min)
bash setup.sh

# 2. submit the 400 jobs: 4 job arrays of 100
bash submit_lsf.sh        # LSF (bsub)
bash submit_slurm.sh      # SLURM (sbatch)
bash run_one_machine.sh   # no scheduler

# 3. progress; resubmit what is missing; collect
bash check_progress.sh
bash submit_lsf.sh --missing       # or submit_slurm.sh --missing
bash collect_results.sh
```

Edit `config.sh` first:

- `R_SETUP` must make both R and a C++ compiler available, e.g. `module load R/4.4.0 gcc/11`.
- `QUEUE`, and `ACCOUNT` if your cluster needs one.

Run `setup.sh` before submitting. It compiles the engine once into `code/cpp/build`; the jobs only load it. If the source file is changed after setup, run `setup.sh` again before submitting. Otherwise every job would compile it at the same time.

`check_progress.sh` and `collect_results.sh` work as in the R version.

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

`NCORES` is used twice: for the bootstrap refits, and for the prediction threads of the C++ engine.

## Time and memory per job

These figures were measured on a Mac (Apple silicon, R 4.6.1) for one replicate at the full settings (10,000 iterations, K = 200, 650 kept draws), with four jobs, each with 4 cores, running at once.

| | R version | C++ version |
|---|---|---|
| Scenario A, Steps 1–3 of the SBART method | 44 min | 26 min |
| Scenario D, Steps 1–3 of the SBART method | 41 min | 30 min |
| Peak memory, A / D | 8.2 / 5.1 GB | 6.5 / 3.7 GB |

On the same Mac with fewer jobs running, using `NCORES=6`, a Scenario A replicate at the full MCMC settings but with 40/16 bootstrap refits (`verify/run_comparison.sh full`) took 34 min in R and 17 min in C++. Of this, Steps 1–3 took 32 min in R and 15 min in C++.

The C++ version is faster in the part after burn-in: 650 times, the posterior survival curve is predicted on the 15 × 512 × 101 AMSE grid, and these predictions run on `NCORES` threads.

The forest updates themselves, about 6 minutes per 3,000 iterations, are the same SoftBart code in both versions and run on one core.

The bootstrap refits are the same R code in both versions. The 500 RSF refits take about 2 hours of core time, divided over `NCORES`.

On the cluster used for the paper, a replicate took much longer than on the Mac: the median time per replicate was 353, 350, 345 and 230 minutes for Scenarios A–D (`../outputs/section4_numbers.txt`). `WALLTIME=12:00` and `MEM_MB=16000` leave a margin.

## Outputs (`collect_results.sh`)

For each scenario, `results/scenX/` gets:

- `RESULTS.txt`: AMSE summary, intervals and importance-sampling ESS;
- `amse_by_rep.csv`;
- `AES_X.pdf` and `AMSE_X.pdf`;
- `intervals_by_rep.rds`.

These come from `code/sim_aggregate.R`, the per-scenario aggregator `scenA_aggregate.R`. It skips the PNG copies of the figures when R cannot write PNG files; the PDFs are always written.

For the paper, `paper_outputs/` gets:

| File | In the paper |
|---|---|
| `Figure2_AMSE_boxplots.pdf` | Figure 2: per-dataset AMSE, one panel per scenario |
| `FigureS1_AES_A.pdf` … `FigureS4_AES_D.pdf` | Figures S1–S4; band = 2.5%–97.5% range of the replicate AES curves |
| `FigureS*_AES_*_ownband.pdf` | not in the paper; band = each method's own 95% interval, averaged over replicates |
| `TableS1_coverage.csv`, `.tex` | Table S1: coverage of each method's own 95% pointwise interval |
| `amse_summary.csv`, `amse_by_replicate.csv` | AMSE mean, SD, median and quartiles |
| `section4_numbers.txt` | the numbers quoted in Section 4, with replicate counts, censoring, settings and run times |

## Verification

| Script | What it checks |
|---|---|
| `verify/check_truncnorm.R` | the engine's truncated-normal sampler gives the same 302,000 draws as `truncnorm::rtruncnorm` |
| `verify/check_forest_identity.R` | the engine's copy of the SoftBart forest gives the same draws as the installed SoftBart package (200 trees, 300 sweeps) |
| `verify/run_comparison.sh quick` | scenario A, replicate 1, short chains: R and C++ versions, every saved field compared |
| `verify/run_comparison.sh small` | scenarios A–D, replicates 1–3, short chains |
| `verify/run_comparison.sh full` | scenario A, replicate 1, the full MCMC settings of `config.sh` |

`verify/sim_worker_R.R` is the R version's `sim_worker.R`, included here for the comparison.

## One detail the engine reproduces

Called from R, SoftBart's `forest$do_gibbs()` does not save R's random-number state, so the next R draw in the iteration starts from the state before the forest update. The engine reproduces this behaviour (with `PutRNGstate()`/`GetRNGstate()` around the forest update) so that its draws match the R version exactly.

## If compiling fails

- **"library 'emutls_w' not found" (a Mac without gfortran):** create a file containing the single line `FLIBS=` and set `R_MAKEVARS_USER` to its path.
- **"undefined reference to pthread_create":** R did not pass `-pthread`. `load_engine.R` adds it to `PKG_LIBS`; check that your `R_SETUP` loads the same compiler that R was built with.
