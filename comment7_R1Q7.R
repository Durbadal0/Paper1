################################################################################
# Comment 7 (Referee 1, Question 7): Sensitivity to CAR Assumption for M
#
# For AA-D:
# (1) Fit spatial CAR model for M_i (our method)
# (2) Fit spatially independent model for M_i (iid Normal prior on logit(M_i))
# (3) For both: run full Step 2 SBART survival analysis
# (4) Compare: Var[logit(M_i)|D_0] vs log(n_{0i})
# (5) Compare: posterior survival / LYG between the two approaches
#
# Usage: Rscript comment7_R1Q7.R
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

cat("=== Comment 7: Spatial vs Independent Model for M_i ===\n\n")

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

m0 <- brfss_data$m
n0 <- brfss_data$n
log_n0 <- log(n0)

# FCR survival data - AA-D stratum
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

cat("AA-D: N =", nrow(surv_data), "subjects,", sum(surv_data$delta), "events\n\n")

################################################################################
# Model 1: Spatial CAR model (our method)
################################################################################

cat("=== Step 1a: Fitting Spatial CAR Model for M ===\n")

fit_spatial <- fit_brfss_car(
  m0 = m0, n0 = n0, A = A,
  num_iter = 5000, burn_in = 2500
)

keep_sp <- (fit_spatial$burn_in + 1):nrow(fit_spatial$M_samples)
M_hat_spatial <- colMeans(fit_spatial$M_samples[keep_sp, ])

# Posterior samples of logit(M_i)
logitM_samples_spatial <- qlogis(fit_spatial$M_samples[keep_sp, ])
# Handle Inf/-Inf from M=0 or M=1
logitM_samples_spatial[!is.finite(logitM_samples_spatial)] <- NA
logitM_var_spatial <- apply(logitM_samples_spatial, 2, var, na.rm = TRUE)

################################################################################
# Model 2: Independent model (no spatial association)
# logit(M_i) ~ iid N(0, sigma0_sq) — same MCMC framework, different prior
################################################################################

cat("\n=== Step 1b: Fitting Independent (iid Normal) Model for M ===\n")

num_iter_indep <- 5000
burn_in_indep <- 2500

M_samples_indep <- matrix(NA, nrow = num_iter_indep, ncol = N)
sigma0_sq_indep <- 1.0
logit_M_indep <- qlogis(pmin(pmax(m0 / n0, 0.01), 0.99))

for (iter in 1:num_iter_indep) {
  # --- Update logit(M) using Metropolis-Hastings ---
  for (i in 1:N) {
    logit_M_prop <- logit_M_indep[i] + rnorm(1, 0, 0.3)

    M_prop_i <- plogis(logit_M_prop)
    M_curr_i <- plogis(logit_M_indep[i])

    # Log likelihood (binomial)
    ll_prop <- dbinom(m0[i], n0[i], M_prop_i, log = TRUE)
    ll_curr <- dbinom(m0[i], n0[i], M_curr_i, log = TRUE)

    # Log prior: iid Normal(0, sigma0_sq) — NOT CAR
    lp_prop <- dnorm(logit_M_prop, 0, sqrt(sigma0_sq_indep), log = TRUE)
    lp_curr <- dnorm(logit_M_indep[i], 0, sqrt(sigma0_sq_indep), log = TRUE)

    log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
    if (log(runif(1)) < log_alpha) {
      logit_M_indep[i] <- logit_M_prop
    }
  }

  # --- Update sigma0_sq (Inverse-Gamma full conditional) ---
  quad_form <- sum(logit_M_indep^2)
  a_post_sig <- 1 + N / 2
  b_post_sig <- 1 + quad_form / 2
  sigma0_sq_indep <- 1 / rgamma(1, a_post_sig, b_post_sig)

  # No rho0 update — no spatial correlation in independent model

  # Store M samples
  M_samples_indep[iter, ] <- plogis(logit_M_indep)

  if (iter %% 1000 == 0) cat("  Independent model: Iteration", iter, "\n")
}

keep_ind <- (burn_in_indep + 1):num_iter_indep
M_hat_indep <- colMeans(M_samples_indep[keep_ind, ])

# Posterior samples of logit(M_i)
logitM_samples_indep <- qlogis(M_samples_indep[keep_ind, ])
logitM_samples_indep[!is.finite(logitM_samples_indep)] <- NA
logitM_var_indep <- apply(logitM_samples_indep, 2, var, na.rm = TRUE)

################################################################################
# Step 2: SBART Survival with Spatial M_hat
################################################################################

cat("\n=== Step 2a: SBART Survival with Spatial CAR M_hat ===\n")
t_start <- Sys.time()

fit_surv_spatial <- sbart_survival_mcmc(
  surv_data = surv_data, A = A, M_hat = M_hat_spatial,
  num_iter = 8000, burn_in = 4000, num_tree = 20
)

t_spatial <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("  Time:", round(t_spatial, 1), "minutes\n")

################################################################################
# Step 2: SBART Survival with Independent M_hat
################################################################################

cat("\n=== Step 2b: SBART Survival with Independent M_hat ===\n")
t_start <- Sys.time()

fit_surv_indep <- sbart_survival_mcmc(
  surv_data = surv_data, A = A, M_hat = M_hat_indep,
  num_iter = 8000, burn_in = 4000, num_tree = 20
)

t_indep <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat("  Time:", round(t_indep, 1), "minutes\n")

################################################################################
# Step 3: Importance Weights
################################################################################

cat("\n=== Step 3a: Importance Weights (Spatial) ===\n")
weights_spatial <- compute_importance_weights(fit_surv_spatial, fit_spatial$M_samples)

cat("\n=== Step 3b: Importance Weights (Independent) ===\n")
weights_indep <- compute_importance_weights(fit_surv_indep, M_samples_indep)

################################################################################
# Compute LYG(10) for TD effect under both models
################################################################################

cat("\n=== Computing LYG(10) comparison ===\n")

compute_lyg_for_model <- function(fit_step2, M_hat_vec, effect = "TD",
                                   tau = 10, n_post = 200, weights = NULL) {
  burn_in <- fit_step2$burn_in
  keep <- (burn_in + 1):length(fit_step2$lambda0_samples)
  max_time <- fit_step2$max_time
  forest <- fit_step2$forest
  cov_names <- fit_step2$cov_names
  X_ref <- fit_step2$X_scaled_ref

  Age_col <- which(cov_names == "Age")
  HR_col <- which(cov_names == "HR")
  Age_med <- median(X_ref[, Age_col])
  HR_med <- median(X_ref[, HR_col])

  # Baseline: median age, HR, Grade3=1, TD=0, BD=0
  # x_covs order: Age, HR, Grade2, Grade3, TD, BD
  x_base <- c(Age_med, HR_med, 0, 1, 0, 0)
  x_alt  <- c(Age_med, HR_med, 0, 1, 0, 0)

  if (effect == "TD") {
    x_alt[5] <- 1   # TD=1
  } else if (effect == "BD") {
    x_alt[6] <- 1   # BD=1
  } else if (effect == "SMU") {
    # For SMU, we vary the M_hat value, not a covariate
    # Will handle separately below
  }

  n_post <- min(n_post, length(keep))
  if (!is.null(weights)) {
    w <- weights[1:length(keep)]
    w[is.na(w)] <- 0
    if (sum(w) > 0) { w <- w / sum(w) } else { w <- rep(1/length(keep), length(keep)) }
    sample_idx <- sample(keep, n_post, replace = TRUE, prob = w)
  } else {
    sample_idx <- sample(keep, n_post)
  }

  # Use average W across counties for a representative patient
  lyg_samples <- numeric(n_post)

  for (s in 1:n_post) {
    iter <- sample_idx[s]
    lam <- fit_step2$lambda0_samples[iter]
    W_s <- fit_step2$W_samples[iter, ]
    W_avg <- median(W_s[W_s > 0 & is.finite(W_s)], na.rm = TRUE)

    if (is.na(lam) || is.na(W_avg) || W_avg <= 0) {
      lyg_samples[s] <- NA; next
    }

    M_avg <- mean(M_hat_vec, na.rm = TRUE)

    n_times <- 40
    times_grid <- seq(0.1, tau, length.out = n_times)
    S_base_vals <- numeric(n_times)
    S_alt_vals <- numeric(n_times)

    for (k in 1:n_times) {
      t_k <- times_grid[k]
      n_int <- 30
      t_int <- seq(0.001, t_k, length.out = n_int)
      dt_int <- t_int[2] - t_int[1]

      if (effect == "SMU") {
        M_q75 <- quantile(M_hat_vec, 0.75, na.rm = TRUE)
        M_q25 <- quantile(M_hat_vec, 0.25, na.rm = TRUE)

        X_base_grid <- cbind(t_int/max_time, rep(M_q75, n_int),
                             matrix(rep(x_base, n_int), nrow=n_int, byrow=TRUE))
        X_alt_grid <- cbind(t_int/max_time, rep(M_q25, n_int),
                            matrix(rep(x_base, n_int), nrow=n_int, byrow=TRUE))
      } else {
        X_base_grid <- cbind(t_int/max_time, rep(M_avg, n_int),
                             matrix(rep(x_base, n_int), nrow=n_int, byrow=TRUE))
        X_alt_grid <- cbind(t_int/max_time, rep(M_avg, n_int),
                            matrix(rep(x_alt, n_int), nrow=n_int, byrow=TRUE))
      }

      b_base <- tryCatch(forest$do_predict(X_base_grid),
                         error = function(e) rep(0, n_int))
      b_alt <- tryCatch(forest$do_predict(X_alt_grid),
                        error = function(e) rep(0, n_int))

      H_base <- lam * W_avg * (sum(pnorm(b_base)) - 0.5*(pnorm(b_base[1]) +
                  pnorm(b_base[n_int]))) * dt_int
      H_alt <- lam * W_avg * (sum(pnorm(b_alt)) - 0.5*(pnorm(b_alt[1]) +
                  pnorm(b_alt[n_int]))) * dt_int

      S_base_vals[k] <- exp(-max(0, H_base))
      S_alt_vals[k] <- exp(-max(0, H_alt))
    }

    # Trapezoidal LYG
    diff_S <- S_base_vals - S_alt_vals
    dt <- diff(times_grid)
    lyg_samples[s] <- sum(dt * (diff_S[-n_times] + diff_S[-1]) / 2)
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

# LYG(10) for TD effect
cat("  LYG(10) for TD effect (spatial)...\n")
lyg_td_spatial <- compute_lyg_for_model(fit_surv_spatial, M_hat_spatial,
                                         effect = "TD", tau = 10, weights = weights_spatial)
cat("  LYG(10) for TD effect (independent)...\n")
lyg_td_indep <- compute_lyg_for_model(fit_surv_indep, M_hat_indep,
                                       effect = "TD", tau = 10, weights = weights_indep)

# LYG(10) for BD effect
cat("  LYG(10) for BD effect (spatial)...\n")
lyg_bd_spatial <- compute_lyg_for_model(fit_surv_spatial, M_hat_spatial,
                                         effect = "BD", tau = 10, weights = weights_spatial)
cat("  LYG(10) for BD effect (independent)...\n")
lyg_bd_indep <- compute_lyg_for_model(fit_surv_indep, M_hat_indep,
                                       effect = "BD", tau = 10, weights = weights_indep)

# LYG(10) for SMU effect
cat("  LYG(10) for SMU effect (spatial)...\n")
lyg_smu_spatial <- compute_lyg_for_model(fit_surv_spatial, M_hat_spatial,
                                          effect = "SMU", tau = 10, weights = weights_spatial)
cat("  LYG(10) for SMU effect (independent)...\n")
lyg_smu_indep <- compute_lyg_for_model(fit_surv_indep, M_hat_indep,
                                        effect = "SMU", tau = 10, weights = weights_indep)

################################################################################
# Compare key parameters between models
################################################################################

cat("\n=== Parameter Comparison ===\n")

keep_sp2 <- (fit_surv_spatial$burn_in + 1):length(fit_surv_spatial$lambda0_samples)
keep_in2 <- (fit_surv_indep$burn_in + 1):length(fit_surv_indep$lambda0_samples)

cat("Spatial model:\n")
cat(sprintf("  lambda_0: %.3f (%.3f, %.3f)\n",
    mean(fit_surv_spatial$lambda0_samples[keep_sp2], na.rm=TRUE),
    quantile(fit_surv_spatial$lambda0_samples[keep_sp2], 0.025, na.rm=TRUE),
    quantile(fit_surv_spatial$lambda0_samples[keep_sp2], 0.975, na.rm=TRUE)))
cat(sprintf("  sigma_1^2: %.3f\n",
    mean(fit_surv_spatial$sigma1_sq_samples[keep_sp2], na.rm=TRUE)))
cat(sprintf("  rho_1: %.3f\n",
    mean(fit_surv_spatial$rho1_samples[keep_sp2], na.rm=TRUE)))

cat("\nIndependent model:\n")
cat(sprintf("  lambda_0: %.3f (%.3f, %.3f)\n",
    mean(fit_surv_indep$lambda0_samples[keep_in2], na.rm=TRUE),
    quantile(fit_surv_indep$lambda0_samples[keep_in2], 0.025, na.rm=TRUE),
    quantile(fit_surv_indep$lambda0_samples[keep_in2], 0.975, na.rm=TRUE)))
cat(sprintf("  sigma_1^2: %.3f\n",
    mean(fit_surv_indep$sigma1_sq_samples[keep_in2], na.rm=TRUE)))
cat(sprintf("  rho_1: %.3f\n",
    mean(fit_surv_indep$rho1_samples[keep_in2], na.rm=TRUE)))

################################################################################
# PLOTS
################################################################################

dir.create("results/comment7", showWarnings = FALSE, recursive = TRUE)

# ---- Plot 1: Var[logit(M_i)|D_0] vs log(n_{0i}) - Spatial CAR only ----
pdf("results/comment7/var_logitM_spatial.pdf", width = 6, height = 5)
par(mar = c(4.5, 5, 0.5, 0.5), family = "serif")

plot(log_n0, logitM_var_spatial, pch = 16, cex = 0.9, col = "black",
     xlab = expression(log(n["0i"])),
     ylab = expression(widehat(Var)(logit(M[i]) ~ "|" ~ D[0])),
     frame.plot = TRUE)

dev.off()
cat("Saved: results/comment7/var_logitM_spatial.pdf\n")

# ---- Plot 2: Comparison - Spatial vs Independent ----
pdf("results/comment7/var_logitM_comparison.pdf", width = 6, height = 5)
par(mar = c(4.5, 5, 0.5, 0.5), family = "serif")

y_range <- range(c(logitM_var_spatial, logitM_var_indep), na.rm = TRUE)

plot(log_n0, logitM_var_spatial, pch = 16, cex = 0.9, col = "black",
     xlab = expression(log(n["0i"])),
     ylab = expression(widehat(Var)(logit(M[i]) ~ "|" ~ D[0])),
     ylim = y_range, frame.plot = TRUE)

points(log_n0, logitM_var_indep, pch = 17, cex = 0.9, col = "firebrick")

dev.off()
cat("Saved: results/comment7/var_logitM_comparison.pdf\n")

# ---- Plot 3: LYG(10) comparison between spatial and independent ----
pdf("results/comment7/lyg10_comparison.pdf", width = 6, height = 5)
par(mar = c(5, 4.5, 0.5, 0.5), family = "serif")

effects <- c("TD", "BD", "SMU")
means_sp <- c(lyg_td_spatial$mean, lyg_bd_spatial$mean, lyg_smu_spatial$mean)
lower_sp <- c(lyg_td_spatial$lower, lyg_bd_spatial$lower, lyg_smu_spatial$lower)
upper_sp <- c(lyg_td_spatial$upper, lyg_bd_spatial$upper, lyg_smu_spatial$upper)

means_in <- c(lyg_td_indep$mean, lyg_bd_indep$mean, lyg_smu_indep$mean)
lower_in <- c(lyg_td_indep$lower, lyg_bd_indep$lower, lyg_smu_indep$lower)
upper_in <- c(lyg_td_indep$upper, lyg_bd_indep$upper, lyg_smu_indep$upper)

y_range <- range(c(lower_sp, upper_sp, lower_in, upper_in), na.rm = TRUE)
if (diff(y_range) < 1e-6) y_range <- y_range + c(-0.1, 0.1)

x_pos <- 1:3
offset <- 0.12

plot(x_pos - offset, means_sp, pch = 16, cex = 1.3, col = "black",
     xlim = c(0.5, 3.5), ylim = y_range,
     xlab = "", ylab = expression(widehat(LYG)(10)),
     xaxt = "n", frame.plot = TRUE)
axis(1, at = x_pos, labels = effects)

arrows(x_pos - offset, lower_sp, x_pos - offset, upper_sp,
       length = 0.04, angle = 90, code = 3, col = "black", lwd = 1.2)

points(x_pos + offset, means_in, pch = 17, cex = 1.3, col = "firebrick")
arrows(x_pos + offset, lower_in, x_pos + offset, upper_in,
       length = 0.04, angle = 90, code = 3, col = "firebrick", lwd = 1.2)

dev.off()
cat("Saved: results/comment7/lyg10_comparison.pdf\n")

# ---- Plot 4: Posterior M_hat comparison ----
pdf("results/comment7/Mhat_comparison.pdf", width = 6, height = 5)
par(mar = c(4.5, 4.5, 0.5, 0.5), family = "serif")

plot(M_hat_spatial, M_hat_indep, pch = 16, cex = 0.9,
     xlab = expression(hat(M)[i] ~ "(Spatial CAR)"),
     ylab = expression(hat(M)[i] ~ "(Independent)"),
     frame.plot = TRUE)
abline(0, 1, lty = 2, col = "gray50")

dev.off()
cat("Saved: results/comment7/Mhat_comparison.pdf\n")

################################################################################
# Summary Statistics
################################################################################

cat("\n=== Variance Comparison Summary ===\n\n")
cat("Var[logit(M_i)|D_0]:\n")
cat("  Spatial CAR:  Mean =", format(mean(logitM_var_spatial, na.rm=TRUE), digits=4),
    " Range = [", format(min(logitM_var_spatial, na.rm=TRUE), digits=4), ",",
    format(max(logitM_var_spatial, na.rm=TRUE), digits=4), "]\n")
cat("  Independent:  Mean =", format(mean(logitM_var_indep, na.rm=TRUE), digits=4),
    " Range = [", format(min(logitM_var_indep, na.rm=TRUE), digits=4), ",",
    format(max(logitM_var_indep, na.rm=TRUE), digits=4), "]\n\n")

cat("Ratio (Independent / Spatial):\n")
ratio <- logitM_var_indep / logitM_var_spatial
ratio <- ratio[is.finite(ratio)]
cat("  Mean ratio:", format(mean(ratio, na.rm=TRUE), digits=3), "\n")
cat("  Spatial model reduces posterior variance by borrowing strength\n\n")

cat("=== LYG(10) Comparison ===\n")
cat(sprintf("TD effect:  Spatial = %.4f (%.4f, %.4f)  |  Independent = %.4f (%.4f, %.4f)\n",
    lyg_td_spatial$mean, lyg_td_spatial$lower, lyg_td_spatial$upper,
    lyg_td_indep$mean, lyg_td_indep$lower, lyg_td_indep$upper))
cat(sprintf("BD effect:  Spatial = %.4f (%.4f, %.4f)  |  Independent = %.4f (%.4f, %.4f)\n",
    lyg_bd_spatial$mean, lyg_bd_spatial$lower, lyg_bd_spatial$upper,
    lyg_bd_indep$mean, lyg_bd_indep$lower, lyg_bd_indep$upper))
cat(sprintf("SMU effect: Spatial = %.4f (%.4f, %.4f)  |  Independent = %.4f (%.4f, %.4f)\n",
    lyg_smu_spatial$mean, lyg_smu_spatial$lower, lyg_smu_spatial$upper,
    lyg_smu_indep$mean, lyg_smu_indep$lower, lyg_smu_indep$upper))

################################################################################
# VIP Comparison
################################################################################

vip_spatial <- compute_vip(fit_surv_spatial)
vip_indep <- compute_vip(fit_surv_indep)

cat("\n=== VIP Comparison ===\n")
cat("Spatial:     ", paste(names(vip_spatial), "=", round(vip_spatial, 3), collapse = ", "), "\n")
cat("Independent: ", paste(names(vip_indep), "=", round(vip_indep, 3), collapse = ", "), "\n")

cat("\n=== Comment 7 Complete ===\n")
