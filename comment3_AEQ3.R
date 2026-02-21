################################################################################
# Comment 3 (AE Question 3): Analysis Excluding Stage Stratification
#
# The AE suggests looking at estimated effects of BD, TD, and SMU from a model
# that excludes the mediator (Stage), to get the total effect of each variable.
#
# Approach:
# (1) Pool all AA subjects (Local + Regional + Distant) = ~6,913 subjects
# (2) Randomly subsample 500 subjects (repeat 2 times with different seeds)
# (3) Fit the 3-step SBART spatial survival model (same as paper, no stage strata)
# (4) Compute LYG(10) for TD, BD, and SMU effects at median covariates
# (5) Compare with the stratified results from the paper
#
# Usage: Rscript comment3_AEQ3.R
################################################################################

rm(list = ls())
set.seed(2024)
setwd("/Users/dghosh35/Desktop/Paper1")

library(SoftBart)
library(MASS)
library(truncnorm)
library(Matrix)
library(readxl)

source("sbart_spatial_survival_correct.R")

cat("=== Comment 3: Analysis Excluding Stage Stratification (AA only) ===\n\n")

################################################################################
# Load Data
################################################################################

A <- as.matrix(read.csv("W.mat.csv", row.names = 1))
N <- nrow(A)

# BRFSS data for AA
brfss_aa_raw <- read_excel("mammo_county_16A.xlsx", skip = 4, col_names = FALSE)
brfss_aa <- data.frame(
  county_raw = as.numeric(brfss_aa_raw[[1]]),
  m0 = as.numeric(brfss_aa_raw[[2]]),
  m1 = as.numeric(brfss_aa_raw[[3]])
)
brfss_aa <- brfss_aa[!is.na(brfss_aa$county_raw), ]
brfss_aa$county <- (brfss_aa$county_raw + 1) / 2
brfss_aa$m <- brfss_aa$m1
brfss_aa$n <- brfss_aa$m0 + brfss_aa$m1

brfss_data <- data.frame(county = 1:N, m = 0, n = 1)
for (i in 1:nrow(brfss_aa)) {
  c_idx <- round(brfss_aa$county[i])
  if (!is.na(c_idx) && c_idx >= 1 && c_idx <= N) {
    brfss_data$m[c_idx] <- brfss_aa$m[i]
    brfss_data$n[c_idx] <- brfss_aa$n[i]
  }
}
brfss_data$n[brfss_data$n == 0] <- 1

# FCR data - ALL AA subjects (Race == 2), pooled across stages
surv_raw <- read.csv("Surv_data2.csv")
names(surv_raw) <- c("id", "time_days", "county", "TX_Delay", "death",
                     "BX_Delay", "Age", "HR_p", "Tgrade", "Race",
                     "SMUprob", "Stage")

aa_all <- surv_raw[surv_raw$Race == 2, ]
cat("Total AA subjects (all stages):", nrow(aa_all), "\n")
cat("  AA-Local:", sum(aa_all$Stage == 1), "\n")
cat("  AA-Regional:", sum(aa_all$Stage == 2), "\n")
cat("  AA-Distant:", sum(aa_all$Stage == 3), "\n\n")

################################################################################
# Step 1: Fit CAR model for M (shared across both subsamples)
################################################################################

cat("=== Step 1: Fitting CAR Model for M ===\n")

fit_step1 <- fit_brfss_car(
  m0 = brfss_data$m, n0 = brfss_data$n, A = A,
  num_iter = 5000, burn_in = 2500
)

keep1 <- (fit_step1$burn_in + 1):nrow(fit_step1$M_samples)
M_hat <- colMeans(fit_step1$M_samples[keep1, ])

################################################################################
# Function to prepare subsample and fit model
################################################################################

prepare_surv_data <- function(raw_data) {
  surv_data <- data.frame(
    county = raw_data$county,
    time = pmax(raw_data$time_days / 365.25, 0.01),
    delta = raw_data$death,
    Age = (raw_data$Age - min(raw_data$Age, na.rm = TRUE)) /
          max(1, max(raw_data$Age, na.rm = TRUE) - min(raw_data$Age, na.rm = TRUE)),
    HR = raw_data$HR_p,
    Grade2 = as.integer(raw_data$Tgrade == 2),
    Grade3 = as.integer(raw_data$Tgrade == 3),
    TD = as.integer(raw_data$TX_Delay == 3),
    BD = as.integer(raw_data$BX_Delay == 3)
  )
  surv_data$time <- pmin(surv_data$time, 15)
  return(surv_data)
}

# Function to compute LYG(10) for a covariate effect
compute_lyg10 <- function(fit_step2, M_hat, x_base, x_alt,
                          n_post_samples = 300) {
  burn_in <- fit_step2$burn_in
  keep <- (burn_in + 1):length(fit_step2$lambda0_samples)
  max_time <- fit_step2$max_time
  forest <- fit_step2$forest

  n_post <- min(n_post_samples, length(keep))
  sample_idx <- sample(keep, n_post)

  # Use county-averaged W (overall effect, not county-specific)
  lyg_samples <- numeric(n_post)

  for (s in 1:n_post) {
    iter <- sample_idx[s]
    lam <- fit_step2$lambda0_samples[iter]
    W_s <- fit_step2$W_samples[iter, ]

    if (is.na(lam)) { lyg_samples[s] <- NA; next }

    # Average W across counties with data
    W_avg <- mean(W_s, na.rm = TRUE)

    # Time grid for integration
    n_grid <- 50
    times_grid <- seq(0.1, 10, length.out = n_grid)
    dt <- times_grid[2] - times_grid[1]

    S_base <- numeric(n_grid)
    S_alt <- numeric(n_grid)

    for (k in 1:n_grid) {
      t_k <- times_grid[k]
      # Integration grid for cumulative hazard
      n_int <- 30
      t_int <- seq(0.001, t_k, length.out = n_int)
      dt_int <- t_int[2] - t_int[1]

      # Baseline profile
      X_base <- cbind(t_int / max_time,
                      rep(x_base[1], n_int),  # M
                      matrix(rep(x_base[-1], n_int), nrow = n_int, byrow = TRUE))
      b_base <- tryCatch(forest$do_predict(X_base),
                         error = function(e) rep(0, n_int))
      H_base <- lam * W_avg * sum(pnorm(b_base)) * dt_int
      S_base[k] <- exp(-H_base)

      # Alternative profile
      X_alt <- cbind(t_int / max_time,
                     rep(x_alt[1], n_int),  # M
                     matrix(rep(x_alt[-1], n_int), nrow = n_int, byrow = TRUE))
      b_alt <- tryCatch(forest$do_predict(X_alt),
                        error = function(e) rep(0, n_int))
      H_alt <- lam * W_avg * sum(pnorm(b_alt)) * dt_int
      S_alt[k] <- exp(-H_alt)
    }

    lyg_samples[s] <- sum(S_base - S_alt) * dt
  }

  lyg_samples <- lyg_samples[!is.na(lyg_samples)]
  return(list(
    mean = mean(lyg_samples),
    median = median(lyg_samples),
    lower = quantile(lyg_samples, 0.025),
    upper = quantile(lyg_samples, 0.975),
    samples = lyg_samples
  ))
}

################################################################################
# Run 2 Subsamples
################################################################################

n_subsample <- 500
n_repeats <- 2
seeds <- c(42, 123)

results_list <- list()

dir.create("results/comment3", showWarnings = FALSE, recursive = TRUE)

for (rep_id in 1:n_repeats) {
  cat("\n================================================================\n")
  cat("=== Subsample", rep_id, "of", n_repeats, "(seed =", seeds[rep_id], ") ===\n")
  cat("================================================================\n\n")

  set.seed(seeds[rep_id])

  # Random subsample of 500 AA subjects (all stages pooled)
  idx <- sample(1:nrow(aa_all), n_subsample)
  sub_data <- aa_all[idx, ]

  cat("Subsample stage distribution:\n")
  cat("  Local:", sum(sub_data$Stage == 1), "\n")
  cat("  Regional:", sum(sub_data$Stage == 2), "\n")
  cat("  Distant:", sum(sub_data$Stage == 3), "\n")
  cat("  Events:", sum(sub_data$death), "\n\n")

  surv_sub <- prepare_surv_data(sub_data)

  # Step 2: SBART survival
  cat("=== Fitting SBART survival model ===\n")
  t_start <- Sys.time()

  fit_step2 <- sbart_survival_mcmc(
    surv_data = surv_sub, A = A, M_hat = M_hat,
    num_iter = 4000, burn_in = 2000, num_tree = 20
  )

  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  cat("  Elapsed:", round(t_elapsed, 1), "minutes\n\n")

  # Compute VIP
  vip <- compute_vip(fit_step2)
  cat("Variable Importance Proportions:\n")
  print(round(vip, 3))
  cat("\n")

  # Median covariate values for the subsample
  Age_med <- median(surv_sub$Age)
  HR_med <- median(surv_sub$HR)
  # Most common grade
  grade_tab <- table(surv_sub$Grade2, surv_sub$Grade3)
  M_med <- mean(M_hat, na.rm = TRUE)

  # Covariate order: M, Age, HR, Grade2, Grade3, TD, BD
  # Baseline profile: median age, HR=1 (positive), Grade3=1, TD=0, BD=0
  x_base_common <- c(M_med, Age_med, HR_med, 0, 1, 0, 0)

  # --- LYG(10) for Treatment Delay ---
  cat("Computing LYG(10) for Treatment Delay effect...\n")
  x_base_td <- x_base_common               # TD=0
  x_alt_td <- x_base_common
  x_alt_td[6] <- 1                         # TD=1
  lyg_td <- compute_lyg10(fit_step2, M_hat, x_base_td, x_alt_td)

  # --- LYG(10) for Biopsy Delay ---
  cat("Computing LYG(10) for Biopsy Delay effect...\n")
  x_base_bd <- x_base_common               # BD=0
  x_alt_bd <- x_base_common
  x_alt_bd[7] <- 1                         # BD=1
  lyg_bd <- compute_lyg10(fit_step2, M_hat, x_base_bd, x_alt_bd)

  # --- LYG(10) for SMU effect ---
  cat("Computing LYG(10) for SMU effect...\n")
  # Compare M at 25th percentile vs 75th percentile
  M_q25 <- quantile(M_hat, 0.25, na.rm = TRUE)
  M_q75 <- quantile(M_hat, 0.75, na.rm = TRUE)
  x_base_smu <- x_base_common
  x_base_smu[1] <- M_q75                   # Higher SMU (better screening)
  x_alt_smu <- x_base_common
  x_alt_smu[1] <- M_q25                    # Lower SMU
  lyg_smu <- compute_lyg10(fit_step2, M_hat, x_base_smu, x_alt_smu)

  results_list[[rep_id]] <- list(
    seed = seeds[rep_id],
    n_local = sum(sub_data$Stage == 1),
    n_regional = sum(sub_data$Stage == 2),
    n_distant = sum(sub_data$Stage == 3),
    n_events = sum(sub_data$death),
    vip = vip,
    lyg_td = lyg_td,
    lyg_bd = lyg_bd,
    lyg_smu = lyg_smu,
    time_min = t_elapsed
  )

  cat("\n--- Results for Subsample", rep_id, "---\n")
  cat(sprintf("  LYG(10) for TD: %.3f (%.3f, %.3f)\n",
              lyg_td$mean, lyg_td$lower, lyg_td$upper))
  cat(sprintf("  LYG(10) for BD: %.3f (%.3f, %.3f)\n",
              lyg_bd$mean, lyg_bd$lower, lyg_bd$upper))
  cat(sprintf("  LYG(10) for SMU (Q75 vs Q25): %.3f (%.3f, %.3f)\n",
              lyg_smu$mean, lyg_smu$lower, lyg_smu$upper))
}

################################################################################
# Summary Comparison
################################################################################

cat("\n\n================================================================\n")
cat("=== SUMMARY: Unstratified AA Model (pooled across stages) ===\n")
cat("================================================================\n\n")

cat("This analysis pools all AA subjects (Local + Regional + Distant)\n")
cat("without conditioning on cancer stage, to estimate the total\n")
cat("effect of TD, BD, and SMU on survival.\n\n")

for (rep_id in 1:n_repeats) {
  r <- results_list[[rep_id]]
  cat(sprintf("Subsample %d (seed=%d): N=%d (L=%d, R=%d, D=%d), events=%d\n",
              rep_id, r$seed, n_subsample,
              r$n_local, r$n_regional, r$n_distant, r$n_events))
  cat(sprintf("  LYG(10) TD:  %.3f (%.3f, %.3f)\n",
              r$lyg_td$mean, r$lyg_td$lower, r$lyg_td$upper))
  cat(sprintf("  LYG(10) BD:  %.3f (%.3f, %.3f)\n",
              r$lyg_bd$mean, r$lyg_bd$lower, r$lyg_bd$upper))
  cat(sprintf("  LYG(10) SMU: %.3f (%.3f, %.3f)\n",
              r$lyg_smu$mean, r$lyg_smu$lower, r$lyg_smu$upper))
  cat(sprintf("  Time: %.1f min\n\n", r$time_min))
}

# Average across repeats
avg_td <- mean(sapply(results_list, function(r) r$lyg_td$mean))
avg_bd <- mean(sapply(results_list, function(r) r$lyg_bd$mean))
avg_smu <- mean(sapply(results_list, function(r) r$lyg_smu$mean))

cat("Average LYG(10) across 2 subsamples:\n")
cat(sprintf("  TD:  %.3f\n", avg_td))
cat(sprintf("  BD:  %.3f\n", avg_bd))
cat(sprintf("  SMU: %.3f\n\n", avg_smu))

################################################################################
# VIP Comparison Plot
################################################################################

# Average VIP across repeats
vip1 <- results_list[[1]]$vip
vip2 <- results_list[[2]]$vip
vip_avg <- (vip1 + vip2) / 2

pdf("results/comment3/vip_unstratified.pdf", width = 6, height = 4.5)
par(mar = c(6, 4.5, 0.5, 0.5), family = "serif")
bp <- barplot(vip_avg, names.arg = names(vip_avg), las = 2,
              ylab = "VIP", col = "gray70", border = "gray30",
              ylim = c(0, max(vip_avg) * 1.2))
abline(h = 1/length(vip_avg), lty = 2, col = "gray50")
dev.off()
cat("Saved: results/comment3/vip_unstratified.pdf\n")

################################################################################
# LYG(10) Comparison Plot: Subsample 1 vs 2
################################################################################

pdf("results/comment3/lyg10_comparison.pdf", width = 6, height = 5)
par(mar = c(5, 4.5, 0.5, 0.5), family = "serif")

effects <- c("TD", "BD", "SMU")
n_eff <- length(effects)

# Collect means and CIs
means_1 <- c(results_list[[1]]$lyg_td$mean,
             results_list[[1]]$lyg_bd$mean,
             results_list[[1]]$lyg_smu$mean)
lower_1 <- c(results_list[[1]]$lyg_td$lower,
             results_list[[1]]$lyg_bd$lower,
             results_list[[1]]$lyg_smu$lower)
upper_1 <- c(results_list[[1]]$lyg_td$upper,
             results_list[[1]]$lyg_bd$upper,
             results_list[[1]]$lyg_smu$upper)

means_2 <- c(results_list[[2]]$lyg_td$mean,
             results_list[[2]]$lyg_bd$mean,
             results_list[[2]]$lyg_smu$mean)
lower_2 <- c(results_list[[2]]$lyg_td$lower,
             results_list[[2]]$lyg_bd$lower,
             results_list[[2]]$lyg_smu$lower)
upper_2 <- c(results_list[[2]]$lyg_td$upper,
             results_list[[2]]$lyg_bd$upper,
             results_list[[2]]$lyg_smu$upper)

y_range <- range(c(lower_1, upper_1, lower_2, upper_2), na.rm = TRUE)

x_pos <- 1:n_eff
offset <- 0.15

plot(x_pos - offset, means_1, pch = 16, cex = 1.3, col = "black",
     xlim = c(0.5, n_eff + 0.5), ylim = y_range,
     xlab = "", ylab = expression(widehat(LYG)(10)),
     xaxt = "n", frame.plot = TRUE)
axis(1, at = x_pos, labels = effects)

arrows(x_pos - offset, lower_1, x_pos - offset, upper_1,
       length = 0.05, angle = 90, code = 3, col = "black", lwd = 1.2)

points(x_pos + offset, means_2, pch = 17, cex = 1.3, col = "gray40")
arrows(x_pos + offset, lower_2, x_pos + offset, upper_2,
       length = 0.05, angle = 90, code = 3, col = "gray40", lwd = 1.2)

dev.off()
cat("Saved: results/comment3/lyg10_comparison.pdf\n")

################################################################################
# Save Results
################################################################################

save(results_list, M_hat, file = "results/comment3/comment3_results.RData")

# Text summary
sink("results/comment3/summary_report.txt")
cat("================================================================\n")
cat("RESPONSE TO AE COMMENT 3\n")
cat("Analysis Excluding Stage Stratification (AA Race Only)\n")
cat(paste("Generated:", Sys.time()), "\n")
cat("================================================================\n\n")

cat("Approach:\n")
cat("  - Pooled all AA subjects (Local + Regional + Distant) into a\n")
cat("    single stratum, removing the stage-based stratification.\n")
cat("  - Randomly subsampled 500 subjects (repeated twice for stability).\n")
cat("  - Fitted the same 3-step SBART spatial survival model as in\n")
cat("    the paper, but without conditioning on cancer stage.\n")
cat("  - Computed LYG(10) for TD, BD, and SMU effects.\n\n")

cat("This estimates the TOTAL effect of each variable on survival,\n")
cat("without conditioning on the potential mediator (Stage).\n\n")

for (rep_id in 1:n_repeats) {
  r <- results_list[[rep_id]]
  cat(sprintf("Subsample %d (seed=%d):\n", rep_id, r$seed))
  cat(sprintf("  Stage composition: L=%d, R=%d, D=%d\n",
              r$n_local, r$n_regional, r$n_distant))
  cat(sprintf("  Events: %d / %d (%.1f%%)\n",
              r$n_events, n_subsample, 100 * r$n_events / n_subsample))
  cat(sprintf("  LYG(10) for TD:  %.3f (%.3f, %.3f)\n",
              r$lyg_td$mean, r$lyg_td$lower, r$lyg_td$upper))
  cat(sprintf("  LYG(10) for BD:  %.3f (%.3f, %.3f)\n",
              r$lyg_bd$mean, r$lyg_bd$lower, r$lyg_bd$upper))
  cat(sprintf("  LYG(10) for SMU: %.3f (%.3f, %.3f)\n",
              r$lyg_smu$mean, r$lyg_smu$lower, r$lyg_smu$upper))
  cat(sprintf("  VIP: %s\n", paste(names(r$vip), "=",
              round(r$vip, 3), collapse = ", ")))
  cat(sprintf("  Computation time: %.1f minutes\n\n", r$time_min))
}

cat("Average across subsamples:\n")
cat(sprintf("  LYG(10) TD:  %.3f\n", avg_td))
cat(sprintf("  LYG(10) BD:  %.3f\n", avg_bd))
cat(sprintf("  LYG(10) SMU: %.3f\n", avg_smu))
cat("\n================================================================\n")
sink()

cat("\nSaved: results/comment3/summary_report.txt\n")
cat("\n=== Comment 3 Complete ===\n")
