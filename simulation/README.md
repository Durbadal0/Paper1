# Simulation study (Section 4 of the main text, Section 1 of the Supplementary Materials)

## Contents

| Folder | Contents |
|---|---|
| `R_version/` | The simulation in plain R: one job per simulated dataset, with scripts for LSF, SLURM or a single machine. |
| `cpp_version/` | The same simulation with the SBART sampler (Steps 1–3) in C++. It uses the same seeds and the same random draws, and gave results identical to the R version in our checks (`cpp_version/VERIFICATION.md`). |
| `outputs/` | The results for the paper: 4 scenarios × 100 datasets, produced with `cpp_version/`. |

## Outputs and the paper

| File in `outputs/` | In the paper |
|---|---|
| `Figure2_AMSE_boxplots.pdf` | Figure 2 of the main text: AMSE of every dataset, one panel per scenario |
| `FigureS1_AES_A.pdf` … `FigureS4_AES_D.pdf` | Figures S1–S4: true survival function and the AES curves of the three methods for county 1. The bands are the 2.5%–97.5% range of the 100 replicate estimates. |
| `TableS1_coverage.csv`, `TableS1_coverage.tex` | Table S1: how often each method's own 95% pointwise interval contains the true survival function, over 100 time points in (0, 3] for the four (x, M1) values of Figures S1–S4, averaged over the 100 datasets |
| `amse_by_replicate.csv` | The AMSE of each method for each of the 400 datasets (the data behind Figure 2) |
| `amse_summary.csv` | Mean, SD, median, quartiles, minimum and maximum of the AMSE by scenario and method |
| `section4_numbers.txt` | Summary per scenario: number of datasets, sample size, censoring, MCMC settings, median run time, AMSE summaries, number of datasets in which SBART has the lowest AMSE, coverage |
| `FigureS*_AES_*_ownband.pdf` | Not in the paper: versions of Figures S1–S4 in which the band is each method's own 95% interval, averaged over the datasets |
| `report.html` | A one-page HTML summary of the same results |
| `build.log` | Log of the cluster job that ran `make_paper_outputs.R` |

All files except `report.html` and `build.log` are written by `cpp_version/code/make_paper_outputs.R` (called by `cpp_version/collect_results.sh`). It builds them from the 400 replicate files `cpp_version/results/scenX/repNN.rds`, which are not included in this repository.

## The design

| | Scenarios A–C | Scenario D |
|---|---|---|
| Counties (clusters) | 15 Florida panhandle counties, real adjacency | 10 |
| Subjects | 1,250 | 1,000 |
| BRFSS-type survey sizes n_0i | 8 to 25 | 8 to 25 |
| t_max for AMSE | 10 | 8 |

The design has the following features:

- **Covariates.** Each subject has p = 15 covariates from Uniform(0, 1). Only x1–x3 enter the true hazard.
- **Censoring.** Censoring times are drawn above the median survival time and calibrated to 12% censoring in every dataset.
- **Scenarios.**
  - A: exponential-link hazard, independent CAR models for the frailties and for M.
  - B: as A, with a time-varying frailty effect.
  - C: as A, with M tied to the frailty: logit(M_i) ~ N(0.9 log W_i, 1 − 0.9²).
  - D: as C with 10 counties and a probit-link hazard.
- **Seeds.** Replicate r uses the seed 20260814 + r.

## Settings used for the paper (`config.sh`)

| Setting | Value |
|---|---|
| Datasets per scenario | 100 |
| SBART MCMC, Step 2 | 10,000 iterations, 3,500 burn-in, every 10th draw kept (650 draws), K = 200 trees |
| Step 1 (CAR model for M) | 12,000 iterations, 6,000 burn-in |
| Algorithm | 3 steps (`STEP3=1`) |
| PH-frailty and RSF intervals | 500 bootstrap refits each |
| RSF | 1,000 trees (250 in each bootstrap refit) |

The outputs were produced on an HPC cluster (LSF, R 4.4.0, GCC 9.5.0) with `cpp_version/`. The median run time per dataset was 353, 350, 345 and 230 minutes for Scenarios A–D.

## How to reproduce

Both versions work the same way. From inside `R_version/` or `cpp_version/`:

```bash
bash setup.sh             # once: R packages, a short test (and, for C++, compiling and checking the engine)
bash submit_lsf.sh        # or submit_slurm.sh, or run_one_machine.sh: 4 job arrays of 100
bash check_progress.sh    # progress; resubmit missing replicates with --missing
bash collect_results.sh   # per-scenario summaries and the paper outputs
```

Edit `config.sh` first: how R is loaded on the cluster (`R_SETUP`), the queue, and the cores and memory per job. The README in each folder gives the details, run times and memory needs.

**Exact numbers depend on the computing environment.** On one machine, the two versions give the same numbers to the last digit. On a different machine, the compiler and its C++ library, the BLAS, and the R and package versions can change the random forests of the RSF and the long MCMC chains, and so the later digits.

For example, for replicate 1 of Scenarios A and D, a run on a Mac and the run on the cluster gave the same PH-frailty AMSE to all printed digits, but slightly different SBART and RSF values.
