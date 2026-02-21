################################################################################
# Comment 4 (Referee 1, Question 4): Variable Cluster Sizes
#
# (1) Plot posterior CI of M_i versus log(n_{0i}) for AA
# (2) Fit SBART survival model with 10,000 iterations
# (3) Plot LYG(10) CI for TD, BD, SMU effects for selected counties
#     vs cluster sizes (AA-D)
# (4) Plot CI width of W_i vs cluster size
#
# Usage: Rscript comment4_R1Q4.R
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

cat("=== Comment 4: Variable Cluster Sizes ===\n\n")

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

# FCR data - AA-D stratum
surv_raw <- read.csv("Surv_data2.csv")
names(surv_raw) <- c("id", "time_days", "county", "TX_Delay", "death",
                     "BX_Delay", "Age", "HR_p", "Tgrade", "Race",
                     "SMUprob", "Stage")

stratum_data <- surv_raw[surv_raw$Race == 2 & surv_raw$Stage == 3, ]
surv_data <- data.frame(
  county = stratum_data$county,
  time = pmax(stratum_data$time_days / 365.25, 0.01),
  delta = stratum_data$death,
  Age = (stratum_data$Age - min(stratum_data$Age)) /
        max(1, max(stratum_data$Age) - min(stratum_data$Age)),
  HR = stratum_data$HR_p,
  Grade2 = as.integer(stratum_data$Tgrade == 2),
  Grade3 = as.integer(stratum_data$Tgrade == 3),
  TD = as.integer(stratum_data$TX_Delay == 3),
  BD = as.integer(stratum_data$BX_Delay == 3)
)
surv_data$time <- pmin(surv_data$time, 15)

cat("AA-D stratum: N =", nrow(surv_data), "\n")
cat("BRFSS counties with data:", sum(brfss_data$n > 1), "\n\n")

################################################################################
# PART 1: Posterior CI of M_i versus log(n_{0i})
################################################################################

cat("=== Part 1: Fitting CAR model for M (Step 1) ===\n")

fit_step1 <- fit_brfss_car(
  m0 = brfss_data$m,
  n0 = brfss_data$n,
  A = A,
  num_iter = 5000,
  burn_in = 2500
)

keep <- (fit_step1$burn_in + 1):nrow(fit_step1$M_samples)
M_hat <- colMeans(fit_step1$M_samples[keep, ])
M_lower <- apply(fit_step1$M_samples[keep, ], 2, quantile, 0.025)
M_upper <- apply(fit_step1$M_samples[keep, ], 2, quantile, 0.975)

n0 <- brfss_data$n
log_n0 <- log(n0)

# Plot 1: Posterior CI of M_i vs log(n_{0i})
dir.create("results/comment4", showWarnings = FALSE, recursive = TRUE)

pdf("results/comment4/Mi_vs_log_n0i.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), family = "serif")

plot(log_n0, M_hat, pch = 16, cex = 0.8,
     xlab = expression(log(n["0i"])), ylab = expression(hat(M)[i]),
     ylim = range(c(M_lower, M_upper)),
     frame.plot = TRUE)

arrows(log_n0, M_lower, log_n0, M_upper,
       length = 0.03, angle = 90, code = 3, col = "gray50", lwd = 0.8)

points(log_n0, M_hat, pch = 16, cex = 0.8, col = "black")

dev.off()
cat("Saved: results/comment4/Mi_vs_log_n0i.pdf\n")

################################################################################
# PART 2: SBART Survival Model (10,000 iterations)
################################################################################

cat("\n=== Part 2: Fitting SBART model (Step 2, 10000 iterations) ===\n")

fit_step2 <- sbart_survival_mcmc(
  surv_data = surv_data,
  A = A,
  M_hat = M_hat,
  num_iter = 10000,
  burn_in = 5000,
  num_tree = 20
)

################################################################################
# Step 3: Importance Weights
################################################################################

cat("\n=== Step 3: Computing Importance Weights ===\n")
weights <- compute_importance_weights(fit_step2, fit_step1$M_samples)

################################################################################
# PART 3: LYG(10) for TD, BD, SMU effects by county
################################################################################

cat("\n=== Part 3: Computing LYG(10) for selected counties ===\n")

# County cluster sizes in AA-D
county_counts <- table(surv_data$county)
county_sizes <- data.frame(
  county = as.integer(names(county_counts)),
  n_i = as.integer(county_counts)
)
county_sizes <- county_sizes[order(county_sizes$n_i), ]

cat("\nCluster size distribution:\n")
print(summary(county_sizes$n_i))

# Select 2-3 small, 1 medium, 2 large
n_cs <- nrow(county_sizes)
small_idx <- county_sizes$county[1:min(3, n_cs)]
med_idx <- county_sizes$county[ceiling(n_cs / 2)]
large_idx <- county_sizes$county[(n_cs - 1):n_cs]
selected_counties <- c(small_idx, med_idx, large_idx)

cat("\nSelected counties:\n")
for (cc in selected_counties) {
  cat("  County", cc, ": n_i =", county_sizes$n_i[county_sizes$county == cc], "\n")
}

# Setup for LYG computation
burn_in <- fit_step2$burn_in
keep_s2 <- (burn_in + 1):length(fit_step2$lambda0_samples)
max_time <- fit_step2$max_time
forest <- fit_step2$forest
cov_names <- fit_step2$cov_names
X_ref <- fit_step2$X_scaled_ref

Age_col <- which(cov_names == "Age")
HR_col <- which(cov_names == "HR")
Age_med <- median(X_ref[, Age_col])
HR_med <- median(X_ref[, HR_col])

# Function to compute S(t) at a single time point
compute_S_at_t <- function(t_k, lam, W_c, M_val, x_covs, forest, max_time) {
  n_int <- 30
  t_int <- seq(0.001, t_k, length.out = n_int)
  dt_int <- t_int[2] - t_int[1]

  X_grid <- cbind(t_int / max_time, rep(M_val, n_int),
                  matrix(rep(x_covs, n_int), nrow = n_int, byrow = TRUE))

  b_grid <- tryCatch(forest$do_predict(X_grid), error = function(e) rep(0, n_int))
  phi_b <- pnorm(b_grid)
  H <- lam * W_c * (sum(phi_b) - 0.5 * (phi_b[1] + phi_b[n_int])) * dt_int
  return(exp(-max(0, H)))
}

# LYG(10) = integral_0^10 [S(t|x_base) - S(t|x_alt)] dt
# For each county, compute LYG(10) for TD, BD, and SMU effects
# x_covs order: Age, HR, Grade2, Grade3, TD, BD

n_post <- min(300, length(keep_s2))
w_lyg <- weights[1:length(keep_s2)]
w_lyg[is.na(w_lyg)] <- 0
if (sum(w_lyg) > 0) { w_lyg <- w_lyg / sum(w_lyg) } else { w_lyg <- rep(1/length(keep_s2), length(keep_s2)) }
sample_idx <- sample(keep_s2, n_post, replace = TRUE, prob = w_lyg)

# Use median W across all posterior samples for stable estimation
W_post_all <- fit_step2$W_samples[keep_s2, ]

tau <- 10
n_times <- 40
times_grid <- seq(0.1, tau, length.out = n_times)

compute_lyg_county <- function(cc, x_base, x_alt, M_val) {
  lyg_samples <- numeric(n_post)

  for (s in 1:n_post) {
    iter <- sample_idx[s]
    lam <- fit_step2$lambda0_samples[iter]
    W_cc <- W_post_all[s, cc]

    if (is.na(lam) || is.na(W_cc) || !is.finite(W_cc) || W_cc <= 0) {
      lyg_samples[s] <- NA; next
    }

    S_base <- sapply(times_grid, function(t) compute_S_at_t(t, lam, W_cc, M_val, x_base, forest, max_time))
    S_alt  <- sapply(times_grid, function(t) compute_S_at_t(t, lam, W_cc, M_val, x_alt, forest, max_time))

    diff_S <- S_base - S_alt
    dt <- diff(times_grid)
    lyg_samples[s] <- sum(dt * (diff_S[-n_times] + diff_S[-1]) / 2)
  }

  lyg_samples <- lyg_samples[!is.na(lyg_samples)]
  if (length(lyg_samples) < 10) return(list(mean=NA, lower=NA, upper=NA))

  return(list(
    mean = mean(lyg_samples),
    lower = quantile(lyg_samples, 0.025),
    upper = quantile(lyg_samples, 0.975)
  ))
}

# Base covariate profile
x_base_common <- c(Age_med, HR_med, 0, 1, 0, 0)  # Grade3=1, TD=0, BD=0

# --- TD effect ---
x_base_td <- x_base_common        # TD=0
x_alt_td <- x_base_common
x_alt_td[5] <- 1                  # TD=1

# --- BD effect ---
x_base_bd <- x_base_common        # BD=0
x_alt_bd <- x_base_common
x_alt_bd[6] <- 1                  # BD=1

# --- SMU effect ---
M_q75 <- quantile(M_hat, 0.75, na.rm = TRUE)
M_q25 <- quantile(M_hat, 0.25, na.rm = TRUE)

lyg_results <- list()

for (cc in selected_counties) {
  cat("\n  County", cc, "(n_i =",
      county_sizes$n_i[county_sizes$county == cc], "):\n")

  M_cc <- M_hat[cc]
  if (is.na(M_cc)) M_cc <- mean(M_hat, na.rm = TRUE)

  cat("    LYG(10) for TD...\n")
  lyg_td <- compute_lyg_county(cc, x_base_td, x_alt_td, M_cc)

  cat("    LYG(10) for BD...\n")
  lyg_bd <- compute_lyg_county(cc, x_base_bd, x_alt_bd, M_cc)

  cat("    LYG(10) for SMU...\n")
  lyg_smu <- compute_lyg_county(cc, x_base_common, x_base_common, M_q75)
  # For SMU we compare high vs low M, with same x
  lyg_smu_samples <- numeric(n_post)
  for (s in 1:n_post) {
    iter <- sample_idx[s]
    lam <- fit_step2$lambda0_samples[iter]
    W_cc_s <- W_post_all[s, cc]
    if (is.na(lam) || is.na(W_cc_s) || !is.finite(W_cc_s) || W_cc_s <= 0) {
      lyg_smu_samples[s] <- NA; next
    }
    S_high <- sapply(times_grid, function(t) compute_S_at_t(t, lam, W_cc_s, M_q75, x_base_common, forest, max_time))
    S_low  <- sapply(times_grid, function(t) compute_S_at_t(t, lam, W_cc_s, M_q25, x_base_common, forest, max_time))
    diff_S <- S_high - S_low
    dt <- diff(times_grid)
    lyg_smu_samples[s] <- sum(dt * (diff_S[-n_times] + diff_S[-1]) / 2)
  }
  lyg_smu_samples <- lyg_smu_samples[!is.na(lyg_smu_samples)]
  if (length(lyg_smu_samples) >= 10) {
    lyg_smu <- list(mean = mean(lyg_smu_samples),
                    lower = quantile(lyg_smu_samples, 0.025),
                    upper = quantile(lyg_smu_samples, 0.975))
  } else {
    lyg_smu <- list(mean = NA, lower = NA, upper = NA)
  }

  lyg_results[[as.character(cc)]] <- list(
    county = cc,
    n_i = county_sizes$n_i[county_sizes$county == cc],
    td = lyg_td, bd = lyg_bd, smu = lyg_smu
  )

  cat(sprintf("    TD:  %.4f (%.4f, %.4f)\n", lyg_td$mean, lyg_td$lower, lyg_td$upper))
  cat(sprintf("    BD:  %.4f (%.4f, %.4f)\n", lyg_bd$mean, lyg_bd$lower, lyg_bd$upper))
  cat(sprintf("    SMU: %.4f (%.4f, %.4f)\n", lyg_smu$mean, lyg_smu$lower, lyg_smu$upper))
}

################################################################################
# Plot: LYG(10) for each effect vs cluster size
################################################################################

lyg_df <- do.call(rbind, lapply(lyg_results, function(r) {
  data.frame(
    county = r$county, n_i = r$n_i,
    td_mean = r$td$mean, td_lower = r$td$lower, td_upper = r$td$upper,
    bd_mean = r$bd$mean, bd_lower = r$bd$lower, bd_upper = r$bd$upper,
    smu_mean = r$smu$mean, smu_lower = r$smu$lower, smu_upper = r$smu$upper
  )
}))
lyg_df <- lyg_df[order(lyg_df$n_i), ]

# --- LYG(10) for TD vs n_i ---
pdf("results/comment4/LYG10_TD_vs_cluster_size.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), family = "serif")
y_rng <- range(c(lyg_df$td_lower, lyg_df$td_upper), na.rm = TRUE)
if (diff(y_rng) < 1e-8) y_rng <- y_rng + c(-0.1, 0.1)
plot(lyg_df$n_i, lyg_df$td_mean, pch = 16, cex = 1.2,
     xlab = expression(n[i]), ylab = expression(widehat(LYG)(10) ~ "for TD"),
     ylim = y_rng, frame.plot = TRUE)
arrows(lyg_df$n_i, lyg_df$td_lower, lyg_df$n_i, lyg_df$td_upper,
       length = 0.05, angle = 90, code = 3, col = "gray40", lwd = 1.2)
points(lyg_df$n_i, lyg_df$td_mean, pch = 16, cex = 1.2)
dev.off()
cat("Saved: results/comment4/LYG10_TD_vs_cluster_size.pdf\n")

# --- LYG(10) for BD vs n_i ---
pdf("results/comment4/LYG10_BD_vs_cluster_size.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), family = "serif")
y_rng <- range(c(lyg_df$bd_lower, lyg_df$bd_upper), na.rm = TRUE)
if (diff(y_rng) < 1e-8) y_rng <- y_rng + c(-0.1, 0.1)
plot(lyg_df$n_i, lyg_df$bd_mean, pch = 16, cex = 1.2,
     xlab = expression(n[i]), ylab = expression(widehat(LYG)(10) ~ "for BD"),
     ylim = y_rng, frame.plot = TRUE)
arrows(lyg_df$n_i, lyg_df$bd_lower, lyg_df$n_i, lyg_df$bd_upper,
       length = 0.05, angle = 90, code = 3, col = "gray40", lwd = 1.2)
points(lyg_df$n_i, lyg_df$bd_mean, pch = 16, cex = 1.2)
dev.off()
cat("Saved: results/comment4/LYG10_BD_vs_cluster_size.pdf\n")

# --- LYG(10) for SMU vs n_i ---
pdf("results/comment4/LYG10_SMU_vs_cluster_size.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), family = "serif")
y_rng <- range(c(lyg_df$smu_lower, lyg_df$smu_upper), na.rm = TRUE)
if (diff(y_rng) < 1e-8) y_rng <- y_rng + c(-0.1, 0.1)
plot(lyg_df$n_i, lyg_df$smu_mean, pch = 16, cex = 1.2,
     xlab = expression(n[i]), ylab = expression(widehat(LYG)(10) ~ "for SMU"),
     ylim = y_rng, frame.plot = TRUE)
arrows(lyg_df$n_i, lyg_df$smu_lower, lyg_df$n_i, lyg_df$smu_upper,
       length = 0.05, angle = 90, code = 3, col = "gray40", lwd = 1.2)
points(lyg_df$n_i, lyg_df$smu_mean, pch = 16, cex = 1.2)
dev.off()
cat("Saved: results/comment4/LYG10_SMU_vs_cluster_size.pdf\n")

################################################################################
# PART 4: CI width of W_i vs cluster size
################################################################################

cat("\n=== Part 4: Posterior W_i CI width vs cluster size ===\n")

# Weighted quantile helper
weighted_quantile <- function(x, w, probs) {
  valid <- is.finite(x) & !is.na(x) & !is.na(w)
  x <- x[valid]; w <- w[valid]
  if (length(x) == 0) return(rep(NA, length(probs)))
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  w <- w / sum(w)
  cum_w <- cumsum(w)
  sapply(probs, function(p) x[which(cum_w >= p)[1]])
}

# Importance weights aligned with post-burn-in samples
w_post <- weights[1:nrow(W_post_all)]
w_post[is.na(w_post)] <- 0
if (sum(w_post) > 0) { w_post <- w_post / sum(w_post) } else { w_post <- rep(1/length(w_post), length(w_post)) }

W_summary <- data.frame(
  county = county_sizes$county,
  n_i = county_sizes$n_i,
  W_median = NA, W_lower = NA, W_upper = NA, W_ci_width = NA
)

for (k in 1:nrow(W_summary)) {
  cc <- W_summary$county[k]
  w_vals <- W_post_all[, cc]
  valid_idx <- is.finite(w_vals) & !is.na(w_vals)
  if (sum(valid_idx) > 10) {
    wq <- weighted_quantile(w_vals, w_post, c(0.5, 0.025, 0.975))
    W_summary$W_median[k] <- wq[1]
    W_summary$W_lower[k] <- wq[2]
    W_summary$W_upper[k] <- wq[3]
    W_summary$W_ci_width[k] <- W_summary$W_upper[k] - W_summary$W_lower[k]
  }
}
W_summary <- W_summary[!is.na(W_summary$W_ci_width), ]

pdf("results/comment4/Wi_CI_vs_cluster_size.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), family = "serif")
plot(log(W_summary$n_i), W_summary$W_ci_width, pch = 16, cex = 0.9,
     xlab = expression(log(n[i])),
     ylab = expression("95% CI width of " * hat(W)[i]),
     frame.plot = TRUE)
dev.off()
cat("Saved: results/comment4/Wi_CI_vs_cluster_size.pdf\n")

################################################################################
# LYG(10) Computation Details Document
################################################################################

sink("results/comment4/LYG10_computation_details.txt")
cat("================================================================\n")
cat("LYG(10) COMPUTATION DETAILS\n")
cat(paste("Generated:", Sys.time()), "\n")
cat("================================================================\n\n")

cat("DEFINITION:\n")
cat("  LYG(a) = integral_0^a [S(t | x_base, M_i, W_i) - S(t | x_alt, M_i, W_i)] dt\n\n")
cat("  where a = 10 years and the survival function is:\n")
cat("  S(t | W_i, M_i; x) = exp(-W_i * lambda_0 * integral_0^t Phi(b(s, M_i, x)) ds)\n\n")
cat("  - b(.) is the SBART function (nonparametric, learned from data)\n")
cat("  - Phi(.) is the standard normal CDF\n")
cat("  - lambda_0 is the baseline hazard rate\n")
cat("  - W_i = exp(R_i) is the spatial frailty for county i\n")
cat("  - M_i is the county-level screening mammography utilization\n\n")

cat("COMPUTATION STEPS:\n")
cat("  1. Draw posterior samples (lambda_0^(s), W_i^(s), b^(s)) from Step 2 MCMC\n")
cat("  2. Compute importance weights omega*_s (Step 3) to correct for plug-in M_hat\n")
cat("  3. For each posterior sample s (drawn with probability proportional to omega*_s):\n")
cat("     a. Set up time grid: t_1,...,t_K from 0.1 to 10 years (K=40)\n")
cat("     b. For each t_k, compute S(t_k) by numerical integration:\n")
cat("        - Set up integration grid: u_1,...,u_J from 0 to t_k (J=30)\n")
cat("        - Evaluate b(u_j, M_i, x) using the SBART forest\n")
cat("        - Compute cumulative hazard H = lambda_0 * W_i * sum[Phi(b(u_j))] * du\n")
cat("          (trapezoidal rule)\n")
cat("        - S(t_k) = exp(-H)\n")
cat("     c. Compute S(t_k) for both baseline (x_base) and alternative (x_alt)\n")
cat("     d. LYG^(s) = integral [S_base - S_alt] dt via trapezoidal rule\n")
cat("  4. Report importance-weighted posterior mean and 95% credible interval of LYG^(s)\n\n")

cat("EFFECTS COMPUTED:\n")
cat("  TD (Treatment Delay):\n")
cat("    x_base: median Age, HR, Grade3=1, TD=0, BD=0\n")
cat("    x_alt:  median Age, HR, Grade3=1, TD=1, BD=0\n")
cat("    LYG(10) = years of life gained if TD is avoided\n\n")

cat("  BD (Biopsy Delay):\n")
cat("    x_base: median Age, HR, Grade3=1, TD=0, BD=0\n")
cat("    x_alt:  median Age, HR, Grade3=1, TD=0, BD=1\n")
cat("    LYG(10) = years of life gained if BD is avoided\n\n")

cat("  SMU (Screening Mammography Utilization):\n")
cat("    x_base and x_alt have same covariates, but M differs:\n")
cat("    S_high uses M at 75th percentile across counties\n")
cat("    S_low  uses M at 25th percentile across counties\n")
cat("    LYG(10) = years of life gained from improving SMU from Q25 to Q75\n\n")

cat("MCMC SETTINGS:\n")
cat("  Step 1 (CAR for M): 5,000 iterations, 2,500 burn-in\n")
cat("  Step 2 (SBART):    10,000 iterations, 5,000 burn-in\n")
cat("  Number of SBART trees: 20\n")
cat("  Posterior samples used for LYG: 300 (randomly sampled from post-burn-in)\n\n")

cat("COUNTY SELECTION:\n")
cat("  Small clusters:  3 counties with smallest n_i in AA-D\n")
cat("  Medium cluster:  1 county at the median n_i\n")
cat("  Large clusters:  2 counties with largest n_i in AA-D\n")
cat("================================================================\n")
sink()

cat("Saved: results/comment4/LYG10_computation_details.txt\n")

# Print full results table
cat("\n=== LYG(10) Results by Cluster Size ===\n")
print(lyg_df)

cat("\n=== Comment 4 Complete ===\n")
