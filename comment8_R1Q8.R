################################################################################
# Comment 8 (Referee 1, Question 8): Computation Time, Hyperparameters,
# Acceptance Rates, Mixing Diagnostics
#
# Uses AA-D stratum. Runs 2 parallel chains with 20,000 iterations each.
#
# Usage: Rscript comment8_R1Q8.R
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

cat("=== Comment 8: Computation, Hyperparameters, Mixing ===\n\n")

################################################################################
# Load Data
################################################################################

A <- as.matrix(read.csv("W.mat.csv", row.names = 1))
N <- nrow(A)

# BRFSS for AA
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

# FCR - AA-D
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
        (max(stratum_data$Age) - min(stratum_data$Age)),
  HR = stratum_data$HR_p,
  Grade2 = as.integer(stratum_data$Tgrade == 2),
  Grade3 = as.integer(stratum_data$Tgrade == 3),
  TD = as.integer(stratum_data$TX_Delay == 3),
  BD = as.integer(stratum_data$BX_Delay == 3)
)
surv_data$time <- pmin(surv_data$time, 15)

n_subjects <- nrow(surv_data)
cat("AA-D: N =", n_subjects, "subjects,", sum(surv_data$delta), "events\n\n")

################################################################################
# MCMC Settings
################################################################################

num_iter_step1 <- 10000
burn_in_step1 <- 5000
num_iter_step2 <- 17000
burn_in_step2 <- 8500
num_tree <- 20
n_chains <- 2

################################################################################
# Step 1 with acceptance tracking (shared across chains)
################################################################################

cat("=== Step 1: CAR for M (with acceptance tracking) ===\n")
step1_start <- Sys.time()

eig <- eigen(A, only.values = TRUE)$values
rho_bounds <- c(1/min(eig), 1/max(eig))

m0 <- brfss_data$m
n0 <- brfss_data$n
M <- (m0 + 0.5) / (n0 + 1)
logit_M <- qlogis(M)
sigma0_sq <- 1
rho0 <- 0

M_samples <- matrix(NA, nrow = num_iter_step1, ncol = N)
sigma0_sq_samples <- rep(NA, num_iter_step1)
rho0_samples <- rep(NA, num_iter_step1)

accept_logitM <- rep(0, N)
accept_rho0 <- 0

for (iter in 1:num_iter_step1) {
  for (i in 1:N) {
    logit_M_prop <- logit_M
    logit_M_prop[i] <- rnorm(1, logit_M[i], 0.3)
    M_prop <- plogis(logit_M_prop)
    M_curr <- plogis(logit_M)
    ll_prop <- dbinom(m0[i], n0[i], M_prop[i], log = TRUE)
    ll_curr <- dbinom(m0[i], n0[i], M_curr[i], log = TRUE)
    lp_prop <- car_log_density(logit_M_prop, A, rho0, sigma0_sq)
    lp_curr <- car_log_density(logit_M, A, rho0, sigma0_sq)
    log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
    if (log(runif(1)) < log_alpha) {
      logit_M <- logit_M_prop
      accept_logitM[i] <- accept_logitM[i] + 1
    }
  }
  M <- plogis(logit_M)

  D <- diag(rowSums(A))
  quad_form <- as.numeric(t(logit_M) %*% (D - rho0 * A) %*% logit_M)
  a_post <- 1 + N/2
  b_post <- 1 + quad_form/2
  sigma0_sq <- 1/rgamma(1, a_post, b_post)

  rho0_prop <- runif(1, max(rho_bounds[1], rho0 - 0.1),
                     min(rho_bounds[2], rho0 + 0.1))
  lp_prop <- car_log_density(logit_M, A, rho0_prop, sigma0_sq)
  lp_curr <- car_log_density(logit_M, A, rho0, sigma0_sq)
  if (log(runif(1)) < (lp_prop - lp_curr)) {
    rho0 <- rho0_prop
    accept_rho0 <- accept_rho0 + 1
  }

  M_samples[iter, ] <- M
  sigma0_sq_samples[iter] <- sigma0_sq
  rho0_samples[iter] <- rho0
  if (iter %% 2000 == 0) cat("  Step 1 iteration", iter, "\n")
}

step1_time <- as.numeric(difftime(Sys.time(), step1_start, units = "mins"))

keep1 <- (burn_in_step1 + 1):num_iter_step1
M_hat <- colMeans(M_samples[keep1, ])

accept_rate_logitM <- mean(accept_logitM / num_iter_step1)
accept_rate_rho0 <- accept_rho0 / num_iter_step1

cat("  Step 1 time:", round(step1_time, 2), "minutes\n")
cat("  Accept rate logit(M_i):", round(accept_rate_logitM * 100, 1), "%\n")
cat("  Accept rate rho_0:", round(accept_rate_rho0 * 100, 1), "%\n\n")

################################################################################
# Step 2: Run 2 chains with different seeds
################################################################################

chain_results <- list()
chain_seeds <- c(42, 123)

for (ch in 1:n_chains) {
  cat("=== Step 2: Chain", ch, "of", n_chains, "(seed =", chain_seeds[ch], ") ===\n")
  set.seed(chain_seeds[ch])

  step2_start <- Sys.time()

  fit_step2 <- sbart_survival_mcmc(
    surv_data = surv_data, A = A, M_hat = M_hat,
    num_iter = num_iter_step2, burn_in = burn_in_step2, num_tree = num_tree
  )

  step2_time <- as.numeric(difftime(Sys.time(), step2_start, units = "mins"))
  cat("  Chain", ch, "time:", round(step2_time, 2), "minutes\n\n")

  chain_results[[ch]] <- list(
    fit = fit_step2,
    time = step2_time
  )

  # Save intermediate results after each chain completes
  save(chain_results, M_hat, step1_time, m0, n0,
       accept_logitM, accept_rho0, sigma0_sq_samples, rho0_samples,
       num_iter_step1, burn_in_step1, num_iter_step2, burn_in_step2,
       num_tree, n_chains, N, n_subjects, accept_rate_logitM, accept_rate_rho0,
       file = "results/comment8/comment8_checkpoint.RData")
  cat("  Saved checkpoint after Chain", ch, "\n")
}

total_time <- step1_time + sum(sapply(chain_results, function(x) x$time))

################################################################################
# Convergence Diagnostics
################################################################################

cat("=== Convergence Diagnostics ===\n\n")

compute_ess <- function(samples) {
  samples <- samples[!is.na(samples)]
  n <- length(samples)
  if (n < 20) return(NA)
  max_lag <- min(n - 1, floor(n / 2), 200)
  acf_result <- acf(samples, lag.max = max_lag, plot = FALSE)
  acf_vals <- acf_result$acf[-1]
  first_neg <- which(acf_vals < 0.05)[1]
  if (is.na(first_neg)) first_neg <- length(acf_vals)
  sum_rho <- sum(acf_vals[1:first_neg])
  tau_int <- 1 + 2 * sum_rho
  return(max(1, n / tau_int))
}

# Extract post-burn-in samples for each chain
keep2 <- (burn_in_step2 + 1):num_iter_step2
lambda0_ch1 <- chain_results[[1]]$fit$lambda0_samples[keep2]
lambda0_ch2 <- chain_results[[2]]$fit$lambda0_samples[keep2]
sigma1_ch1 <- chain_results[[1]]$fit$sigma1_sq_samples[keep2]
sigma1_ch2 <- chain_results[[2]]$fit$sigma1_sq_samples[keep2]
rho1_ch1 <- chain_results[[1]]$fit$rho1_samples[keep2]
rho1_ch2 <- chain_results[[2]]$fit$rho1_samples[keep2]

# Gelman-Rubin R-hat (simplified)
compute_rhat <- function(chain1, chain2) {
  chain1 <- chain1[!is.na(chain1)]
  chain2 <- chain2[!is.na(chain2)]
  n <- min(length(chain1), length(chain2))
  if (n < 20) return(NA)
  chain1 <- chain1[1:n]
  chain2 <- chain2[1:n]

  W <- (var(chain1) + var(chain2)) / 2
  grand_mean <- (mean(chain1) + mean(chain2)) / 2
  B <- n * ((mean(chain1) - grand_mean)^2 + (mean(chain2) - grand_mean)^2)
  V_hat <- (1 - 1/n) * W + (1/n) * B
  return(sqrt(V_hat / W))
}

rhat_lam <- compute_rhat(lambda0_ch1, lambda0_ch2)
rhat_sig <- compute_rhat(sigma1_ch1, sigma1_ch2)
rhat_rho <- compute_rhat(rho1_ch1, rho1_ch2)

# ESS per chain
ess_lam_ch1 <- compute_ess(lambda0_ch1)
ess_lam_ch2 <- compute_ess(lambda0_ch2)
ess_sig_ch1 <- compute_ess(sigma1_ch1)
ess_sig_ch2 <- compute_ess(sigma1_ch2)
ess_rho_ch1 <- compute_ess(rho1_ch1)
ess_rho_ch2 <- compute_ess(rho1_ch2)

cat("R-hat:\n")
cat("  lambda_0:", round(rhat_lam, 3), "\n")
cat("  sigma_1^2:", round(rhat_sig, 3), "\n")
cat("  rho_1:", round(rhat_rho, 3), "\n\n")

cat("ESS (Chain 1 / Chain 2):\n")
cat("  lambda_0:", round(ess_lam_ch1), "/", round(ess_lam_ch2), "\n")
cat("  sigma_1^2:", round(ess_sig_ch1), "/", round(ess_sig_ch2), "\n")
cat("  rho_1:", round(ess_rho_ch1), "/", round(ess_rho_ch2), "\n\n")

################################################################################
# Trace Plots (both chains overlaid)
################################################################################

dir.create("results/comment8", showWarnings = FALSE, recursive = TRUE)

# Helper function for trace plots
plot_trace <- function(filename, samples_ch1, samples_ch2, ylab_expr,
                       burn_in, width = 7, height = 3.5) {
  pdf(filename, width = width, height = height)
  par(mar = c(4, 4.5, 0.5, 0.5), family = "serif")
  n_iter <- length(samples_ch1)
  y_rng <- range(c(samples_ch1, samples_ch2), na.rm = TRUE)
  plot(1:n_iter, samples_ch1, type = "l", col = "gray40",
       xlab = "Iteration", ylab = ylab_expr,
       ylim = y_rng, frame.plot = TRUE)
  lines(1:n_iter, samples_ch2, col = "firebrick", lty = 1)
  abline(v = burn_in, lty = 2, col = "black")
  dev.off()
}

# Step 2 trace plots
plot_trace("results/comment8/trace_lambda0.pdf",
           chain_results[[1]]$fit$lambda0_samples,
           chain_results[[2]]$fit$lambda0_samples,
           expression(lambda[0]), burn_in_step2)

plot_trace("results/comment8/trace_sigma1sq.pdf",
           chain_results[[1]]$fit$sigma1_sq_samples,
           chain_results[[2]]$fit$sigma1_sq_samples,
           expression(sigma[1]^2), burn_in_step2)

plot_trace("results/comment8/trace_rho1.pdf",
           chain_results[[1]]$fit$rho1_samples,
           chain_results[[2]]$fit$rho1_samples,
           expression(rho[1]), burn_in_step2)

cat("Saved trace plots to results/comment8/\n")

# Step 1 trace plots
pdf("results/comment8/trace_sigma0sq.pdf", width = 7, height = 3.5)
par(mar = c(4, 4.5, 0.5, 0.5), family = "serif")
plot(1:num_iter_step1, sigma0_sq_samples, type = "l", col = "gray40",
     xlab = "Iteration", ylab = expression(sigma[0]^2), frame.plot = TRUE)
abline(v = burn_in_step1, lty = 2, col = "black")
dev.off()

pdf("results/comment8/trace_rho0.pdf", width = 7, height = 3.5)
par(mar = c(4, 4.5, 0.5, 0.5), family = "serif")
plot(1:num_iter_step1, rho0_samples, type = "l", col = "gray40",
     xlab = "Iteration", ylab = expression(rho[0]), frame.plot = TRUE)
abline(v = burn_in_step1, lty = 2, col = "black")
dev.off()

################################################################################
# ACF Plots (using chain 1 post-burn-in)
################################################################################

pdf("results/comment8/acf_lambda0.pdf", width = 5, height = 3.5)
par(mar = c(4, 4.5, 0.5, 0.5), family = "serif")
acf(lambda0_ch1[!is.na(lambda0_ch1)], main = "", lag.max = 80,
    ylab = expression("ACF of " * lambda[0]))
dev.off()

pdf("results/comment8/acf_sigma1sq.pdf", width = 5, height = 3.5)
par(mar = c(4, 4.5, 0.5, 0.5), family = "serif")
acf(sigma1_ch1[!is.na(sigma1_ch1)], main = "", lag.max = 80,
    ylab = expression("ACF of " * sigma[1]^2))
dev.off()

pdf("results/comment8/acf_rho1.pdf", width = 5, height = 3.5)
par(mar = c(4, 4.5, 0.5, 0.5), family = "serif")
acf(rho1_ch1[!is.na(rho1_ch1)], main = "", lag.max = 80,
    ylab = expression("ACF of " * rho[1]))
dev.off()

cat("Saved ACF plots\n\n")

################################################################################
# Generate Report
################################################################################

cat("=== REVIEWER REPORT FOR COMMENT 8 ===\n\n")

full_dataset_size <- 76164
avg_step2_time <- mean(sapply(chain_results, function(x) x$time))
time_per_1000 <- (step1_time + avg_step2_time) / (n_subjects / 1000)

report_lines <- c(
  "================================================================",
  "RESPONSE TO REVIEWER 1, COMMENT 8",
  paste("Stratum: AA-D, N =", n_subjects, "subjects"),
  paste("Number of chains:", n_chains),
  paste("Generated:", Sys.time()),
  "================================================================",
  "",
  "1. COMPUTATION TIME",
  paste("   Step 1 (CAR for M,", num_iter_step1, "iter):", round(step1_time, 1), "minutes"),
  paste("   Step 2 Chain 1 (SBART,", num_iter_step2, "iter):",
        round(chain_results[[1]]$time, 1), "minutes"),
  paste("   Step 2 Chain 2 (SBART,", num_iter_step2, "iter):",
        round(chain_results[[2]]$time, 1), "minutes"),
  paste("   Total (1 chain):", round(step1_time + avg_step2_time, 1), "minutes"),
  paste("   Time per 1000 subjects:", round(time_per_1000, 1), "min"),
  paste("   Extrapolated for full dataset (N=76,164):",
        round(time_per_1000 * full_dataset_size / 1000, 0), "min",
        "(~", round(time_per_1000 * full_dataset_size / 1000 / 60, 1), "hours)"),
  "",
  "2. HYPERPARAMETERS",
  "   Step 1 (CAR for M):",
  "     sigma_0^2 ~ Inverse-Gamma(1, 1)",
  "     rho_0 ~ Uniform(1/e_min, 1/e_max)  [eigenvalue-based bounds]",
  "     MH proposal for logit(M_i): N(current, 0.3^2)",
  "     MH proposal for rho_0: Uniform(rho_0 +/- 0.1)",
  "",
  "   Step 2 (SBART survival):",
  "     lambda_0 ~ Gamma(1, 1)",
  "     sigma_1^2 ~ Inverse-Gamma(1, 1)",
  "     rho_1 ~ Uniform(1/e_min, 1/e_max)",
  "     MH proposal for R_i: N(current, 0.2^2)",
  "     MH proposal for rho_1: Uniform(rho_1 +/- 0.05)",
  "",
  "   SBART (SoftBart package):",
  paste("     Number of trees K =", num_tree),
  "     Tree prior: alpha = 0.95, beta = 2",
  paste("     Leaf prior: sigma_mu = 3/(2*sqrt(K)) =", round(3/(2*sqrt(num_tree)), 4)),
  "",
  "3. ACCEPTANCE RATES (MH steps)",
  paste("   logit(M_i):", round(accept_rate_logitM * 100, 1), "%"),
  paste("   rho_0:", round(accept_rate_rho0 * 100, 1), "%"),
  "",
  "4. HYPERPARAMETER TUNING",
  "   MH step sizes were tuned via pilot runs to achieve",
  "   acceptance rates in the 20-50% range.",
  "   SBART hyperparameters use defaults from the SoftBart R package",
  "   (Linero & Yang, 2018), which have been empirically validated.",
  "   CAR priors are weakly informative IG(1,1) to let data dominate.",
  "",
  "5. MIXING DIAGNOSTICS (2 chains)",
  "",
  "   R-hat (Gelman-Rubin):",
  paste("     lambda_0:", round(rhat_lam, 3)),
  paste("     sigma_1^2:", round(rhat_sig, 3)),
  paste("     rho_1:", round(rhat_rho, 3)),
  "",
  "   ESS (Chain 1 / Chain 2):",
  paste("     lambda_0:", round(ess_lam_ch1), "/", round(ess_lam_ch2)),
  paste("     sigma_1^2:", round(ess_sig_ch1), "/", round(ess_sig_ch2)),
  paste("     rho_1:", round(ess_rho_ch1), "/", round(ess_rho_ch2)),
  "",
  "   Posterior means (Chain 1 / Chain 2):",
  paste("     lambda_0:", round(mean(lambda0_ch1, na.rm=TRUE), 3), "/",
        round(mean(lambda0_ch2, na.rm=TRUE), 3)),
  paste("     sigma_1^2:", round(mean(sigma1_ch1, na.rm=TRUE), 3), "/",
        round(mean(sigma1_ch2, na.rm=TRUE), 3)),
  paste("     rho_1:", round(mean(rho1_ch1, na.rm=TRUE), 3), "/",
        round(mean(rho1_ch2, na.rm=TRUE), 3)),
  "",
  "   Trace plots and ACF plots saved in results/comment8/",
  "================================================================"
)

writeLines(report_lines, "results/comment8/reviewer_report.txt")
cat(paste(report_lines, collapse = "\n"))

cat("\n\n=== Comment 8 Complete ===\n")
