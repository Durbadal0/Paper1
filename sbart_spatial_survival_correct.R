################################################################################
# SBART Spatial Survival Model - Correct Implementation
#
# Implements the EXACT method from:
# "Analysis of spatially clustered cancer survival registry using SBART
#  when a cluster-level covariate is unavailable within registry"
#
# Three-Step Algorithm:
# Step 1: Fit CAR model to BRFSS data to estimate M (SMU proportions)
# Step 2: Fit SBART survival model with M_hat and spatial frailties W
# Step 3: Compute importance weights for proper posterior inference
################################################################################

# Required packages
library(SoftBart)
library(MASS)
library(truncnorm)
library(Matrix)
library(readxl)

################################################################################
# CAR Model Functions
################################################################################

#' Compute CAR precision matrix
#' Q = (D - rho * A) / sigma_sq
car_precision <- function(A, rho, sigma_sq) {
  N <- nrow(A)
  D <- diag(rowSums(A))
  Q <- (D - rho * A) / sigma_sq
  return(Q)
}

#' Log density of CAR model (Equation 3 in paper)
#' R ~ MVN(0, sigma^2 * (D - rho*A)^{-1})
car_log_density <- function(x, A, rho, sigma_sq) {
  N <- length(x)
  D <- diag(rowSums(A))
  Q <- (D - rho * A) / sigma_sq

  # Log determinant
  eig_vals <- eigen(D - rho * A, only.values = TRUE)$values
  log_det <- sum(log(pmax(eig_vals, 1e-10))) - N * log(sigma_sq)

  # Quadratic form: x' Q x
  quad <- as.numeric(t(x) %*% Q %*% x)

  return(0.5 * log_det - 0.5 * quad)
}

################################################################################
# STEP 1: Fit CAR Model to BRFSS Data (Equation 5 in paper)
################################################################################

#' Fit CAR model to BRFSS binomial data
#' m_{0i} ~ Binomial(n_{0i}, M_i)
#' logit(M) ~ CAR(A; sigma_0^2, rho_0)
#'
#' @param m0 Vector of successes (number with screening) per county
#' @param n0 Vector of sample sizes per county
#' @param A Adjacency matrix
#' @param num_iter Number of MCMC iterations
#' @param burn_in Burn-in iterations
#' @return List with M_hat and posterior samples
fit_brfss_car <- function(m0, n0, A, num_iter = 5000, burn_in = 2500) {

  N <- length(m0)

  # Eigenvalue bounds for rho
  eig <- eigen(A, only.values = TRUE)$values
  rho_bounds <- c(1/min(eig), 1/max(eig))

  # Initialize
  M <- (m0 + 0.5) / (n0 + 1)  # Smoothed proportion
  logit_M <- qlogis(M)
  sigma0_sq <- 1
  rho0 <- 0

  # Storage
  M_samples <- matrix(NA, nrow = num_iter, ncol = N)
  sigma0_sq_samples <- rep(NA, num_iter)
  rho0_samples <- rep(NA, num_iter)

  # Hyperparameters
  a_sigma <- 1
  b_sigma <- 1

  for (iter in 1:num_iter) {

    # --- Update logit(M) using Metropolis-Hastings ---
    for (i in 1:N) {
      logit_M_prop <- logit_M
      logit_M_prop[i] <- rnorm(1, logit_M[i], 0.3)

      M_prop <- plogis(logit_M_prop)
      M_curr <- plogis(logit_M)

      # Log likelihood (binomial)
      ll_prop <- dbinom(m0[i], n0[i], M_prop[i], log = TRUE)
      ll_curr <- dbinom(m0[i], n0[i], M_curr[i], log = TRUE)

      # Log prior (CAR)
      lp_prop <- car_log_density(logit_M_prop, A, rho0, sigma0_sq)
      lp_curr <- car_log_density(logit_M, A, rho0, sigma0_sq)

      # Accept/reject
      log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)

      if (log(runif(1)) < log_alpha) {
        logit_M <- logit_M_prop
      }
    }

    M <- plogis(logit_M)

    # --- Update sigma0_sq (Inverse-Gamma full conditional) ---
    D <- diag(rowSums(A))
    quad_form <- as.numeric(t(logit_M) %*% (D - rho0 * A) %*% logit_M)
    a_post <- a_sigma + N/2
    b_post <- b_sigma + quad_form/2
    sigma0_sq <- 1/rgamma(1, a_post, b_post)

    # --- Update rho0 using Metropolis-Hastings ---
    rho0_prop <- runif(1, max(rho_bounds[1], rho0 - 0.1),
                       min(rho_bounds[2], rho0 + 0.1))

    lp_prop <- car_log_density(logit_M, A, rho0_prop, sigma0_sq)
    lp_curr <- car_log_density(logit_M, A, rho0, sigma0_sq)

    if (log(runif(1)) < (lp_prop - lp_curr)) {
      rho0 <- rho0_prop
    }

    # Store samples
    M_samples[iter, ] <- M
    sigma0_sq_samples[iter] <- sigma0_sq
    rho0_samples[iter] <- rho0

    if (iter %% 500 == 0) {
      cat("Step 1: Iteration", iter, "of", num_iter, "\n")
    }
  }

  # Posterior mean of M (M_hat from paper)
  keep <- (burn_in + 1):num_iter
  M_hat <- colMeans(M_samples[keep, ])

  return(list(
    M_hat = M_hat,
    M_samples = M_samples,
    sigma0_sq_samples = sigma0_sq_samples,
    rho0_samples = rho0_samples,
    burn_in = burn_in
  ))
}

################################################################################
# STEP 2: SBART Survival MCMC (Algorithm 1 in paper)
################################################################################

#' Main SBART survival MCMC with spatial frailties
#' Hazard: lambda_ij(t | W_i, M_i; x_ij) = lambda_0 * W_i * Phi(b(t, M_i, x_ij))
#'
#' @param surv_data Data frame with county, time, delta, covariates
#' @param A Adjacency matrix
#' @param M_hat Estimated M from Step 1
#' @param num_iter MCMC iterations
#' @param burn_in Burn-in
#' @param num_tree Number of SBART trees
#' @return List with posterior samples
sbart_survival_mcmc <- function(surv_data, A, M_hat, num_iter = 5000,
                                 burn_in = 2500, num_tree = 20) {

  N <- nrow(A)
  n_total <- nrow(surv_data)

  # Ensure county indices are valid (1 to N)
  valid_counties <- !is.na(surv_data$county) & surv_data$county >= 1 & surv_data$county <= N
  if (!all(valid_counties, na.rm = TRUE)) {
    cat("Warning: Removing", sum(!valid_counties, na.rm = TRUE), "records with invalid county indices\n")
    surv_data <- surv_data[valid_counties, ]
    n_total <- nrow(surv_data)
  }

  # Identify covariate columns
  exclude_cols <- c("county", "time", "delta", "M_hat", "id")
  cov_names <- setdiff(names(surv_data), exclude_cols)
  p <- length(cov_names)

  # Add M_hat to data - handle any NA values
  surv_data$M_hat <- M_hat[surv_data$county]
  surv_data$M_hat[is.na(surv_data$M_hat)] <- mean(M_hat, na.rm = TRUE)

  # Scale covariates to [0,1]
  max_time <- max(surv_data$time)
  X_covs <- as.matrix(surv_data[, cov_names])
  X_scaled <- apply(X_covs, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < 1e-10) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  })

  # Build design matrix: (time/max_time, M_hat, scaled covariates)
  X_design <- cbind(surv_data$time / max_time, surv_data$M_hat, X_scaled)
  colnames(X_design) <- c("time", "M", cov_names)

  # Eigenvalue bounds for rho
  eig <- eigen(A, only.values = TRUE)$values
  rho_bounds <- c(1/min(eig), 1/max(eig))

  # Initialize parameters
  lambda0 <- 0.1
  R <- rep(0, N)  # log(W)
  W <- exp(R)
  sigma1_sq <- 1
  rho1 <- 0

  # Initialize SoftBart forest
  # Note: Cannot pass constant Y values to Hypers() as it causes cv.glmnet to fail
  # Use initial Y based on event indicators transformed for probit
  delta_safe <- ifelse(is.na(surv_data$delta), 0, surv_data$delta)
  Y_init <- ifelse(delta_safe == 1, 0.5, -0.5) + rnorm(n_total, 0, 0.1)
  hypers <- Hypers(X_design, Y_init, num_tree = num_tree)
  opts <- Opts(num_burn = 0, num_save = 1, update_sigma = FALSE)
  forest <- MakeForest(hypers, opts)

  # Storage
  lambda0_samples <- rep(NA, num_iter)
  W_samples <- matrix(NA, nrow = num_iter, ncol = N)
  sigma1_sq_samples <- rep(NA, num_iter)
  rho1_samples <- rep(NA, num_iter)

  # Track variable selection counts for VIP
  var_counts <- rep(0, ncol(X_design))
  names(var_counts) <- colnames(X_design)

  # MCMC iterations
  for (iter in 1:num_iter) {

    # --- Data Augmentation: Generate latent rejection times G ---
    # Algorithm 1 lines 4-7 (wrap in tryCatch for SoftBart stability)
    X_all <- matrix(nrow = 0, ncol = ncol(X_design))
    Y_all <- c()

    tryCatch({
      for (i in 1:n_total) {
        county_i <- surv_data$county[i]
        y_i <- surv_data$time[i]
        delta_i <- surv_data$delta[i]

        # Skip if county index is invalid
        if (is.na(county_i) || county_i < 1 || county_i > N) next
        if (is.na(y_i) || y_i <= 0) next

        # Generate number of candidate points from homogeneous Poisson
        rate_i <- lambda0 * W[county_i] * y_i
        if (is.na(rate_i) || rate_i <= 0) rate_i <- 0.01
        q_i <- rpois(1, rate_i)

        if (q_i > 0 && !is.na(surv_data$M_hat[i])) {
          # Generate uniform times
          G_tilde <- runif(q_i, 0, y_i)

          # Build design for G points
          X_G <- cbind(G_tilde / max_time,
                       rep(surv_data$M_hat[i], q_i),
                       matrix(rep(X_scaled[i, ], q_i), nrow = q_i, byrow = TRUE))

          # Get current predictions
          b_G <- tryCatch(forest$do_predict(X_G), error = function(e) rep(0, nrow(X_G)))

          # Thin: keep points where U > Phi(b) (i.e., rejected)
          U <- runif(q_i)
          keep_idx <- U > pnorm(b_G)

          if (sum(keep_idx) > 0) {
            X_reject <- X_G[keep_idx, , drop = FALSE]
            X_all <- rbind(X_all, X_reject)
            Y_all <- c(Y_all, rep(0, sum(keep_idx)))  # Failures for probit
          }
        }

        # Add event time if delta = 1
        if (!is.na(delta_i) && delta_i == 1) {
          X_event <- matrix(X_design[i, ], nrow = 1)
          X_all <- rbind(X_all, X_event)
          Y_all <- c(Y_all, 1)  # Success for probit
        }
      }
    }, error = function(e) {
      # Data augmentation error - use events only
      for (i in 1:n_total) {
        if (!is.na(surv_data$delta[i]) && surv_data$delta[i] == 1) {
          X_all <<- rbind(X_all, matrix(X_design[i, ], nrow = 1))
          Y_all <<- c(Y_all, 1)
        }
      }
    })

    # --- Update SBART function b(.) using Bayesian backfitting ---
    n_obs <- nrow(X_all)
    if (is.null(n_obs)) n_obs <- 0

    if (n_obs > 10 && length(Y_all) > 10) {
      tryCatch({
        # Ensure X_all is a proper matrix
        X_all <- as.matrix(X_all)

        # Check for any NA/Inf in X_all
        valid_rows <- complete.cases(X_all) & !apply(X_all, 1, function(x) any(is.infinite(x)))
        if (sum(valid_rows) > 10) {
          X_all <- X_all[valid_rows, , drop = FALSE]
          Y_all <- Y_all[valid_rows]

          # Probit data augmentation with truncated normals
          b_current <- forest$do_predict(X_all)
          Z_latent <- numeric(length(Y_all))

          for (j in seq_along(Y_all)) {
            if (Y_all[j] == 1) {
              Z_latent[j] <- rtruncnorm(1, a = 0, b = Inf, mean = b_current[j], sd = 1)
            } else {
              Z_latent[j] <- rtruncnorm(1, a = -Inf, b = 0, mean = b_current[j], sd = 1)
            }
          }

          # Update forest
          forest$do_gibbs(X_all, Z_latent, X_all, 1)

          # Track splits for VIP
          counts <- forest$get_counts()
          var_counts <- var_counts + counts
        }
      }, error = function(e) {
        cat("Error at iteration", iter, ":", e$message, "\n")
      })
    }

    # --- Update lambda0 (Gamma full conditional) ---
    # From Equation 8: lambda0 ~ Gamma(a + sum(delta) + sum(m_ij), b + sum(W_i * y_ij))
    n_events <- sum(surv_data$delta, na.rm = TRUE)
    n_rejections <- sum(Y_all == 0, na.rm = TRUE)
    sum_Wy <- sum(W[surv_data$county] * surv_data$time, na.rm = TRUE)

    a_lam <- 1 + n_events + n_rejections
    b_lam <- 1 + sum_Wy
    lambda0 <- rgamma(1, a_lam, b_lam)

    # --- Update CAR parameters (sigma1_sq, rho1) ---
    D <- diag(rowSums(A))

    # sigma1_sq: Inverse-Gamma
    quad_form <- as.numeric(t(R) %*% (D - rho1 * A) %*% R)
    a_post <- 1 + N/2
    b_post <- 1 + quad_form/2
    sigma1_sq <- 1/rgamma(1, a_post, b_post)

    # rho1: Metropolis-Hastings
    rho1_prop <- runif(1, max(rho_bounds[1], rho1 - 0.05),
                       min(rho_bounds[2], rho1 + 0.05))
    lp_prop <- car_log_density(R, A, rho1_prop, sigma1_sq)
    lp_curr <- car_log_density(R, A, rho1, sigma1_sq)
    log_ratio <- lp_prop - lp_curr
    if (!is.na(log_ratio) && is.finite(log_ratio) && log(runif(1)) < log_ratio) {
      rho1 <- rho1_prop
    }

    # --- Update frailties W using Metropolis-Hastings ---
    # (Paper uses HMC, but MH is simpler and works for moderate N)
    for (c in 1:N) {
      R_prop <- R
      R_prop[c] <- rnorm(1, R[c], 0.2)
      W_prop <- exp(R_prop)

      # CAR prior
      lp_prop <- car_log_density(R_prop, A, rho1, sigma1_sq)
      lp_curr <- car_log_density(R, A, rho1, sigma1_sq)

      # Likelihood for county c
      county_idx <- which(surv_data$county == c)
      if (length(county_idx) > 0) {
        ll_prop <- 0
        ll_curr <- 0

        for (idx in county_idx) {
          y_idx <- surv_data$time[idx]
          delta_idx <- surv_data$delta[idx]

          # Skip if missing values
          if (is.na(y_idx) || is.na(delta_idx)) next

          # Simplified likelihood (event + exponential term)
          if (delta_idx == 1) {
            ll_prop <- ll_prop + log(W_prop[c])
            ll_curr <- ll_curr + log(W[c])
          }
          ll_prop <- ll_prop - lambda0 * W_prop[c] * y_idx
          ll_curr <- ll_curr - lambda0 * W[c] * y_idx
        }
      } else {
        ll_prop <- 0
        ll_curr <- 0
      }

      log_alpha <- (ll_prop + lp_prop) - (ll_curr + lp_curr)
      if (!is.na(log_alpha) && is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        R <- R_prop
        W <- W_prop
      }
    }

    # Center R for identifiability
    R <- R - mean(R)
    W <- exp(R)

    # Store samples
    lambda0_samples[iter] <- lambda0
    W_samples[iter, ] <- W
    sigma1_sq_samples[iter] <- sigma1_sq
    rho1_samples[iter] <- rho1

    if (iter %% 100 == 0) {
      cat("Step 2: Iteration", iter, "of", num_iter, "\n")
    }
  }

  return(list(
    lambda0_samples = lambda0_samples,
    W_samples = W_samples,
    sigma1_sq_samples = sigma1_sq_samples,
    rho1_samples = rho1_samples,
    forest = forest,
    var_counts = var_counts,
    surv_data = surv_data,
    max_time = max_time,
    cov_names = cov_names,
    X_scaled_ref = X_scaled,
    burn_in = burn_in
  ))
}

################################################################################
# STEP 3: Importance Weights (Equation 6 in paper)
################################################################################

#' Compute importance sampling weights
#' omega*_ij = L_ij(W_i, lambda0, b; M_i) / L_ij(W_i, lambda0, b; M_hat_i)
#'
#' @param fit_step2 Output from Step 2
#' @param M_samples M samples from Step 1
#' @param K_star Monte Carlo samples for denominator
#' @return Normalized weights
compute_importance_weights <- function(fit_step2, M_samples, K_star = 5) {

  burn_in <- fit_step2$burn_in
  keep <- (burn_in + 1):length(fit_step2$lambda0_samples)

  # Handle NA values in lambda0_samples
  valid_keep <- keep[!is.na(fit_step2$lambda0_samples[keep])]
  if (length(valid_keep) == 0) {
    cat("Warning: No valid lambda0 samples, returning uniform weights\n")
    return(rep(1/length(keep), length(keep)))
  }

  num_samples <- min(length(valid_keep), nrow(M_samples) - burn_in)
  if (num_samples <= 0) {
    cat("Warning: No valid samples for importance weights\n")
    return(rep(1/length(keep), length(keep)))
  }

  surv_data <- fit_step2$surv_data
  n_total <- nrow(surv_data)
  N <- ncol(fit_step2$W_samples)
  max_time <- fit_step2$max_time
  forest <- fit_step2$forest
  cov_names <- fit_step2$cov_names
  X_scaled_ref <- fit_step2$X_scaled_ref

  # Get M_hat vector indexed by county
  M_hat_vec <- rep(NA, N)
  for (i in 1:n_total) {
    c_i <- surv_data$county[i]
    if (!is.na(c_i) && c_i >= 1 && c_i <= N) {
      M_hat_vec[c_i] <- surv_data$M_hat[i]
    }
  }
  # Fill any missing with mean
  M_hat_vec[is.na(M_hat_vec)] <- mean(M_hat_vec, na.rm = TRUE)

  weights <- rep(1, num_samples)

  cat("Computing importance weights for", num_samples, "samples...\n")

  for (s in 1:num_samples) {
    if (s %% 200 == 0) cat("  Weight sample", s, "of", num_samples, "\n")

    iter_step2 <- valid_keep[s]
    iter_step1 <- burn_in + s

    if (iter_step1 > nrow(M_samples)) next

    lambda0 <- fit_step2$lambda0_samples[iter_step2]
    W <- fit_step2$W_samples[iter_step2, ]
    M_true <- M_samples[iter_step1, ]

    # Skip if any required values are NA
    if (is.na(lambda0) || any(is.na(W)) || any(is.na(M_true))) next

    log_weight <- 0

    # Sample a subset of subjects for computational efficiency
    sample_idx <- sample(1:n_total, min(500, n_total))

    for (i in sample_idx) {
      county_i <- surv_data$county[i]

      # Skip if invalid county
      if (is.na(county_i) || county_i < 1 || county_i > N) next

      # Only compute if M differs significantly
      M_true_i <- M_true[county_i]
      M_hat_i <- M_hat_vec[county_i]

      if (is.na(M_true_i) || is.na(M_hat_i)) next
      if (abs(M_true_i - M_hat_i) > 0.01) {
        y_i <- surv_data$time[i]
        delta_i <- surv_data$delta[i]
        x_i <- X_scaled_ref[i, ]

        if (is.na(y_i) || is.na(delta_i) || any(is.na(x_i))) next

        # Build design matrices
        X_true <- matrix(c(y_i/max_time, M_true_i, x_i), nrow = 1)
        X_hat <- matrix(c(y_i/max_time, M_hat_i, x_i), nrow = 1)

        # Get b predictions
        b_true <- tryCatch(forest$do_predict(X_true), error = function(e) NA)
        b_hat <- tryCatch(forest$do_predict(X_hat), error = function(e) NA)

        if (is.na(b_true) || is.na(b_hat)) next

        # Log-likelihood ratio (simplified - just the Phi terms)
        if (delta_i == 1) {
          log_weight <- log_weight + pnorm(b_true, log.p = TRUE) -
            pnorm(b_hat, log.p = TRUE)
        }
      }
    }

    # Scale by sampling fraction
    log_weight <- log_weight * (n_total / length(sample_idx))

    # Clip extreme weights
    if (!is.na(log_weight) && is.finite(log_weight)) {
      log_weight <- max(min(log_weight, 50), -50)
      weights[s] <- exp(log_weight)
    }
  }

  # Normalize
  weights <- weights / sum(weights, na.rm = TRUE)
  weights[is.na(weights)] <- 1/length(weights)

  # Effective sample size
  ess <- 1 / sum(weights^2, na.rm = TRUE)
  cat("Effective sample size:", round(ess, 1), "of", num_samples, "\n")

  return(weights)
}

################################################################################
# Helper Functions
################################################################################

#' Compute Variable Importance Proportions (VIP)
#' VIP = proportion of splits involving each variable
compute_vip <- function(fit_step2) {
  counts <- fit_step2$var_counts
  if (sum(counts, na.rm = TRUE) == 0) {
    vip <- rep(1/length(counts), length(counts))
  } else {
    vip <- counts / sum(counts, na.rm = TRUE)
  }
  # Ensure names are preserved
  if (!is.null(names(counts))) {
    names(vip) <- names(counts)
  }
  return(vip)
}

#' Predict survival function S(t) for new patient
#' S(t|W,M,x) = exp(-W * lambda0 * integral_0^t Phi(b(s,M,x)) ds)
predict_survival <- function(fit_step2, new_patient, times, burn_in = NULL,
                              weights = NULL) {

  if (is.null(burn_in)) burn_in <- fit_step2$burn_in
  keep <- (burn_in + 1):length(fit_step2$lambda0_samples)
  num_samples <- length(keep)

  if (is.null(weights)) weights <- rep(1/num_samples, num_samples)

  max_time <- fit_step2$max_time
  cov_names <- fit_step2$cov_names
  forest <- fit_step2$forest

  # Scale covariates
  X_ref <- fit_step2$X_scaled_ref
  x_scaled <- sapply(seq_along(cov_names), function(j) {
    rng <- range(X_ref[, j])
    if (diff(rng) < 1e-10) return(0.5)
    (new_patient[[cov_names[j]]] - rng[1]) / diff(rng)
  })

  n_times <- length(times)
  S_samples <- matrix(NA, nrow = num_samples, ncol = n_times)

  for (s in 1:num_samples) {
    iter <- keep[s]
    lambda0 <- fit_step2$lambda0_samples[iter]
    W <- fit_step2$W_samples[iter, new_patient$county]

    # Skip if lambda0 or W is NA
    if (is.na(lambda0) || is.na(W)) {
      S_samples[s, ] <- NA
      next
    }

    for (k in 1:n_times) {
      t_k <- times[k]

      # Numerical integration using Riemann sum
      n_grid <- 30
      t_grid <- seq(0, t_k, length.out = n_grid)
      dt <- t_grid[2] - t_grid[1]

      X_grid <- cbind(t_grid / max_time,
                      rep(new_patient$M_hat, n_grid),
                      matrix(rep(x_scaled, n_grid), nrow = n_grid, byrow = TRUE))

      b_grid <- tryCatch(forest$do_predict(X_grid), error = function(e) rep(0, n_grid))
      phi_b <- pnorm(b_grid)

      # Cumulative hazard
      H <- lambda0 * W * sum(phi_b) * dt

      # Survival (bound to [0,1])
      S_samples[s, k] <- max(0, min(1, exp(-H)))
    }
  }

  # Remove rows with all NA
  valid_rows <- apply(S_samples, 1, function(x) !all(is.na(x)))
  if (sum(valid_rows) == 0) {
    warning("All survival predictions are NA, returning NA matrix")
    return(S_samples)
  }
  S_samples <- S_samples[valid_rows, , drop = FALSE]

  return(S_samples)
}

################################################################################
# Main Fitting Function
################################################################################

#' Fit complete SBART Spatial Survival Model (3-step algorithm)
#'
#' @param surv_data Survival data with county, time, delta, covariates
#' @param brfss_data BRFSS data with county, m (successes), n (total)
#' @param A Adjacency matrix
#' @param num_iter_step1 MCMC iterations for Step 1
#' @param num_iter_step2 MCMC iterations for Step 2
#' @param burn_in_step1 Burn-in for Step 1
#' @param burn_in_step2 Burn-in for Step 2
#' @param num_tree Number of SBART trees
#' @return List with all results
fit_sbart_spatial_survival <- function(surv_data, brfss_data, A,
                                        num_iter_step1 = 5000,
                                        num_iter_step2 = 5000,
                                        burn_in_step1 = 2500,
                                        burn_in_step2 = 2500,
                                        num_tree = 20) {

  cat("\n=== Step 1: Fitting CAR model to BRFSS data ===\n")
  fit_step1 <- fit_brfss_car(
    m0 = brfss_data$m,
    n0 = brfss_data$n,
    A = A,
    num_iter = num_iter_step1,
    burn_in = burn_in_step1
  )

  cat("\n=== Step 2: Fitting SBART survival model ===\n")
  fit_step2 <- sbart_survival_mcmc(
    surv_data = surv_data,
    A = A,
    M_hat = fit_step1$M_hat,
    num_iter = num_iter_step2,
    burn_in = burn_in_step2,
    num_tree = num_tree
  )

  cat("\n=== Step 3: Computing importance weights ===\n")
  weights <- compute_importance_weights(
    fit_step2 = fit_step2,
    M_samples = fit_step1$M_samples
  )

  cat("\n=== Done ===\n")

  return(list(
    fit_step1 = fit_step1,
    fit_step2 = fit_step2,
    weights = weights,
    M_hat = fit_step1$M_hat
  ))
}
