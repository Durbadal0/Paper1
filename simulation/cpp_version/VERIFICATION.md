# Verification of the C++ version against the R version

We ran both versions on the same Mac and compared every quantity each one saves in `repNN.rds`. The fields compared include:

- AMSE (whole range and t95 range);
- the AES curves;
- every 95% band, for SBART, PH-frailty and RSF;
- the importance weights and ESS;
- M_hat and W_hat;
- the posterior means of λ0, σ1² and ρ1;
- acceptance rates;
- the traces.

Only the run time was left out. "Identical" means R's `identical()` returns TRUE for the field: the same numbers to the last bit.

The machine was macOS on Apple silicon, with:

- R 4.6.1, using R's reference BLAS/LAPACK;
- SoftBart 1.0.3, truncnorm 1.0.9, ranger 0.18.0, survival 3.8.6, MASS 7.3.65;
- Rcpp 1.1.2 and RcppArmadillo 15.4.2.1.

## Results

| Check | Result |
|---|---|
| Truncated-normal sampler vs `truncnorm::rtruncnorm` (`verify/check_truncnorm.R`) | 302,000 draws on each side: identical |
| Engine's forest code vs installed SoftBart (`verify/check_forest_identity.R`): 200 trees, 17 covariates, 300 sweeps | predictions, split counts and the next random number: identical |
| Scenarios A–D × replicates 1–3; 3-step algorithm; 600 iterations, 300 burn-in, K = 50 (`verify/run_comparison.sh small`) | 12 of 12 replicates identical in every field (`verify/local_results/small_4scenarios_x3reps.txt`) |
| Scenarios B and D, replicate 4; 2-step algorithm (`STEP3=0`) | identical in every field (`verify/local_results/two_step_B_D.txt`) |
| Scenario A, replicate 1; full settings of `config.sh` (10,000 iterations, 3,500 burn-in, K = 200, every 10th kept, NQ = 50); 40/16 bootstrap refits; final files (`verify/run_comparison.sh full`) | identical in every field (`verify/local_results/full_settings_A_D.txt`) |
| Scenarios A and D, replicate 1; the same settings; the earlier engine | identical in every field except the σ1² trace and the posterior mean of σ1² (`s1_post`), which differed by at most 3.6 × 10⁻¹⁵ and 2.2 × 10⁻¹⁶ (see the note below) |
| The R version vs the single-core scripts `scenA.R` and `scenC.R` | identical in every field |

## Notes on the full-settings run

The σ1² difference in the earlier run came from that engine's own loop for the quadratic form t(R) Q R, which R's BLAS rounds slightly differently. The engine now takes this number from R itself.

The following were obtained after that change, and are identical in every field, traces included:

- the 12-replicate comparison (final files);
- the full-settings run (final files);
- the 2-step check, run with the final engine just before a last edit of the R driver script. That edit, which frees memory before the bootstrap, is the same in both versions.

## Checking on another machine

`setup.sh` repeats three of these checks on the machine where it is run: the truncated-normal check, the forest check, and one short replicate run by both versions.

On a different machine, compiler or BLAS, the two versions should still agree. To confirm across all four scenarios, run:

```bash
bash verify/run_comparison.sh small
```
