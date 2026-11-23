#!/usr/bin/env Rscript
################################################################################
# Data analysis of "Analysis of spatially clustered cancer survival registry data
# using SBART when a cluster-level covariate is unavailable" (Section 5 and
# Supplement Section 2): Florida Cancer Registry (FCR) breast cancer, 67 counties,
# 6 strata = race (AA; WA = non-AA) x stage at diagnosis (L, R, D).

## ---- helpers shared by fit and report
## weighted HPD: the shortest interval carrying `level` of the weight
whpd <- function(v, w, level = 0.95) {
  o <- order(v); v <- v[o]; ww <- w[o]/sum(w); cw <- cumsum(ww); n <- length(v)
  j <- 1L; best <- c(v[1], v[n]); bw <- Inf
  for (i in seq_len(n)) {
    lo <- if (i == 1L) 0 else cw[i-1L]
    if (j < i) j <- i
    while (j <= n && cw[j] - lo < level - 1e-12) j <- j + 1L
    if (j > n) break
    if (v[j] - v[i] < bw) { bw <- v[j] - v[i]; best <- c(v[i], v[j]) }
  }
  best }
## weighted quantiles (p may be a vector)
wq <- function(v, w, p) {
  o <- order(v)
  cw <- cumsum(w[o])/sum(w)
  sapply(p, function(pp) v[o][which(cw >= pp)[1]]) }
## FCR extract: the 12 columns of Surv_data2.csv, by position
read_fcr <- function() {
  d <- read.csv("Surv_data2.csv")
  names(d) <- c("id", "days", "county", "TX", "death", "BX", "Age", "HR", "Tg", "Race", "SMUprob", "Stage")
  d }
## BRFSS D_0 for one race: m_0i of n_0i report mammography; a county absent from the file
## gets n_0i = 0 (no binomial term; the CAR prior fills it in)
read_brfss <- function(file, N) {
  bf <- suppressMessages(readxl::read_excel(file, skip = 4, col_names = FALSE))
  bcode <- suppressWarnings(as.numeric(bf[[1]]))
  ## county code = (FIPS + 1)/2; Miami-Dade is FIPS 086 (formerly 025, i.e. county 13)
  bcty <- ifelse(!is.na(bcode) & bcode == 86, 13, (bcode + 1)/2)
  b_no  <- suppressWarnings(as.numeric(bf[[2]]))     # column "0": no mammography
  b_yes <- suppressWarnings(as.numeric(bf[[3]]))     # column "1": mammography
  ok <- !is.na(bcty) & bcty >= 1 & bcty <= N
  m_0 <- n_0i <- rep(0, N)
  m_0[round(bcty[ok])]  <- b_yes[ok]
  n_0i[round(bcty[ok])] <- b_no[ok] + b_yes[ok]
  list(m_0 = m_0, n_0i = n_0i) }
## deviance residuals from the Cox-Snell residuals cs and the event indicators delta
dev_res <- function(delta, cs) {
  mart <- delta - cs
  sign(mart) * sqrt(pmax(-2*(mart + delta*log(pmax(cs, 1e-12))), 0)) }
## VIP: weighted mean share of the splitting rules that use each variable (grade = 2 indicators)
vip_of <- function(vd, w) {
  pr <- vd/pmax(rowSums(vd), 1)
  v <- colSums(w*pr)/sum(w)
  c(Age = v[["Age"]], HR = v[["HR"]], Grade = v[["Grade2"]] + v[["Grade3"]], SMU = v[["M"]],
    TD = v[["TD"]], BD = v[["BD"]], time = v[["time"]]) }

################################################################################
# PART A -- fit: data processing and Steps 1-3 for one stratum
# Writes fcr_<S>.rds (summaries) and fcr_<S>_details.rds (draws pooled by report).
################################################################################
fit_stratum <- function(args) {
  options(warn = -1)
  race   <- toupper(args[1])
  stage  <- toupper(args[2])
  outdir <- args[3]
  argn <- function(k, d) if (length(args) >= k && nzchar(args[k])) as.numeric(args[k]) else d
  num_iter  <- argn(4, 30000); burn_in <- argn(5, 15000); num_tree <- argn(6, 200)
  use_step3 <- argn(7, 1);     THIN    <- argn(8, 10);    NQ       <- argn(9, 50)
  n1        <- argn(10, 30000); b1     <- argn(11, 10000)
  spatial   <- argn(12, 1)          # 1 = CAR frailty; 0 = independent frailty
  county_oi <- argn(13, 6)          # Broward (county codes follow FIPS order)
  t_hor     <- argn(14, 10)
  seed      <- argn(15, 20260816)   # change only to run an additional chain
  stopifnot(race %in% c("AA","WA"), stage %in% c("L","R","D"),
            "numeric arguments must not be NA" = !anyNA(c(num_iter, burn_in, num_tree, use_step3, THIN,
                                                          NQ, n1, b1, spatial, county_oi, t_hor, seed)),
            "use_step3 must be 1 (3-step) or 0 (2-step)" = use_step3 %in% c(0, 1),
            "spatial must be 1 or 0" = spatial %in% c(0, 1),
            THIN >= 1, NQ >= 1, num_tree >= 1, burn_in >= 0, b1 >= 0, b1 < n1 - 1,
            "integer arguments" = all(c(num_iter, burn_in, num_tree, THIN, NQ, n1, b1, county_oi, seed) ==
                                      round(c(num_iter, burn_in, num_tree, THIN, NQ, n1, b1, county_oi, seed))),
            num_iter >= 1, seed >= 1, seed <= .Machine$integer.max,
            "t_horizon must be >= 10 with t = 5 and t = 10 on its grid (e.g. 10, 20, 25, 50)" =
              t_hor >= 10 && abs(500/t_hor - round(500/t_hor)) < 1e-9,
            "fewer than 2 retained draws" = sum(seq_len(num_iter) > burn_in & seq_len(num_iter) %% THIN == 0) >= 2)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  suppressMessages({library(SoftBart); library(truncnorm); library(readxl)})
  logf <- file.path(outdir, "analysis.log")
  lg <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), " [", race, "-", stage, "] ", ...)
                        cat(m, "\n", file = logf, append = TRUE); cat(m, "\n"); flush.console() }
  t_start <- Sys.time()
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20260816)   # Step 1 always uses this stream, so every chain shares the same M_hat

  ## ==== 0. Data
  A <- as.matrix(read.csv("W.mat.csv", row.names = 1))
  N <- nrow(A); D <- diag(rowSums(A))
  stopifnot(county_oi %in% seq_len(N))
  ## check that the adjacency follows the county order of the data (known neighbours)
  adj_ok <- N == 67 && all(A[cbind(c(6, 13, 17, 29, 37), c(50, 44, 57, 52, 20))] == 1)

  surv_raw <- read_fcr()
  race_code  <- if (race == "AA") 2 else 1          # Race: 2 = African American
  SHORT_DELAY <- 1    # TX and BX (treatment and biopsy delay) are coded 1 = short, 3 = long
  stage_code <- switch(stage, L = 1, R = 2, D = 3)
  st <- surv_raw[surv_raw$Race == race_code & surv_raw$Stage == stage_code, ]

  ## follow-up in years
  yr <- st$days / 365.25
  delta <- as.integer(st$death == 1)
  yr <- pmax(yr, 1/365.25)                           # a same-day event is 1 day, not 0

  sd_df <- data.frame(
    county = st$county, time = yr, delta = delta,
    Age    = st$Age,
    HR     = as.integer(st$HR == 1),
    Grade2 = as.integer(st$Tg == 2),
    Grade3 = as.integer(st$Tg == 3),
    TD     = as.integer(st$TX == SHORT_DELAY),   # indicator of a short treatment delay
    BD     = as.integer(st$BX == SHORT_DELAY))   # indicator of a short biopsy delay
  sd_df <- sd_df[!is.na(sd_df$county) & sd_df$county >= 1 & sd_df$county <= N &
                 !is.na(sd_df$time) & sd_df$time > 0 & !is.na(sd_df$delta), ]
  nt <- nrow(sd_df); cty <- sd_df$county

  cov_names <- c("Age","HR","Grade2","Grade3","TD","BD")
  Xraw <- as.matrix(sd_df[, cov_names])
  ## scale continuous covariates to (0,1) as Section 2.2 requires; binary stay 0/1
  X <- apply(Xraw, 2, function(x) { r <- range(x, na.rm = TRUE)
    if (diff(r) < 1e-10) rep(0.5, length(x)) else (x - r[1])/diff(r) })

  ## BRFSS D_0 for this race (m_0i, n_0i; see read_brfss)
  bf_file <- if (race == "AA") "mammo_county_16A.xlsx" else "mammo_county_16W.xlsx"
  bf <- read_brfss(bf_file, N)
  m_0 <- bf$m_0
  n_0i <- bf$n_0i
  has0 <- n_0i > 0
  lg(sprintf("stratum n=%d events=%d (%.1f%% censored) | counties with BRFSS data: %d/%d",
             nt, sum(sd_df$delta), 100*mean(sd_df$delta == 0), sum(has0), N))
  if (!adj_ok) lg("WARNING: W.mat.csv does not match the county order of the data (known Florida ",
                  "neighbours are not adjacent); the CAR models will use wrong neighbours")
  lg(sprintf("settings: K=%d iter=%d burn=%d THIN=%d NQ=%d step3=%d spatial=%d",
             num_tree, num_iter, burn_in, THIN, NQ, use_step3, spatial))

  ## support of rho from the eigenvalues of D^{-1/2} A D^{-1/2}, capped at (-1, 1)
  dv <- rowSums(A); An <- diag(1/sqrt(dv)) %*% A %*% diag(1/sqrt(dv))
  ev_n <- eigen((An + t(An))/2, only.values = TRUE)$values
  rho_bounds <- c(max(-1, 1/min(ev_n)), min(1, 1/max(ev_n)))

  ## log|D - rho A| from the same eigenvalues
  car_ld <- function(x, rho, s2)
    0.5*(sum(log(dv)) + sum(log(pmax(1 - rho*ev_n, 1e-300))) - N*log(s2)) -
    0.5*as.numeric(t(x) %*% (D - rho*A) %*% x)/s2
  ## truncated uniform random walk for rho: proposal and window width (Hastings factor)
  rho_prop <- function(r, h) runif(1, max(rho_bounds[1], r - h), min(rho_bounds[2], r + h))
  rho_win  <- function(r, h) min(rho_bounds[2], r + h) - max(rho_bounds[1], r - h)
  iid_ld <- function(x, rho, s2) sum(dnorm(x, 0, sqrt(s2), log = TRUE))
  frail_ld <- if (spatial == 1) car_ld else iid_ld     # spatial = 0 is the competitor

  ## ==== 1. STEP 1 -- CAR model for M from D_0 alone
  lg("=== STEP 1: CAR for M from BRFSS ===")
  logitM <- qlogis((m_0 + 0.5)/(n_0i + 1)); s0 <- 1; rho0 <- 0.9
  Msave <- matrix(NA_real_, n1, N); acc1 <- rep(0, N); acc_r0 <- 0
  rho0_tr <- s0_tr <- numeric(n1); H0 <- 0.02          # window half-width for rho_0
  for (it in 1:n1) {
    for (i in 1:N) {
      prop <- logitM; prop[i] <- rnorm(1, logitM[i], 0.3)
      ll_p <- if (has0[i]) dbinom(m_0[i], n_0i[i], plogis(prop[i]),   log = TRUE) else 0
      ll_c <- if (has0[i]) dbinom(m_0[i], n_0i[i], plogis(logitM[i]), log = TRUE) else 0
      la <- (car_ld(prop, rho0, s0) + ll_p) - (car_ld(logitM, rho0, s0) + ll_c)
      if (is.finite(la) && log(runif(1)) < la) { logitM <- prop; acc1[i] <- acc1[i] + 1 }
    }
    qf <- as.numeric(t(logitM) %*% (D - rho0*A) %*% logitM)
    s0 <- 1/rgamma(1, 1 + N/2, 1 + qf/2)
    rp <- rho_prop(rho0, H0)
    if (log(runif(1)) < car_ld(logitM, rp, s0) - car_ld(logitM, rho0, s0) +
                        log(rho_win(rho0, H0)) - log(rho_win(rp, H0))) {
      rho0 <- rp; acc_r0 <- acc_r0 + 1 }
    Msave[it, ] <- plogis(logitM); rho0_tr[it] <- rho0; s0_tr[it] <- s0
  }
  M_post <- Msave[(b1+1):n1, , drop = FALSE]
  M_hat  <- colMeans(M_post)
  var_logitM <- apply(qlogis(pmin(pmax(M_post, 1e-6), 1-1e-6)), 2, var)   # Figure S13
  lg(sprintf("Step 1 done | acc logitM=%.3f acc rho0=%.3f | M_hat range %.3f-%.3f",
             mean(acc1)/n1, acc_r0/n1, min(M_hat), max(M_hat)))
  if (seed != 20260816) set.seed(seed)   # an additional chain: new stream for Steps 2-3 only

  ## ==== 2. STEP 2 -- SBART survival with Poisson-thinning data augmentation
  lg("=== STEP 2: SBART survival ===")
  max_time <- max(sd_df$time)
  Xd <- cbind(sd_df$time/max_time, M_hat[cty], X)
  colnames(Xd) <- c("time", "M", cov_names)
  Yinit <- ifelse(sd_df$delta == 1, 0.5, -0.5) + rnorm(nt, 0, 0.1)
  ## SoftBart defaults (Linero and Yang, 2018); sigma = 1 for the probit augmentation,
  ## leaf scale 3/(2 sqrt(K))
  hyp <- Hypers(Xd, Yinit, num_tree = num_tree, sigma_hat = 1, k = 1/3)
  next_stream <- function(st) assign(".Random.seed", parallel::nextRNGStream(st), envir = globalenv())
  st <- .Random.seed
  forest <- MakeForest(hyp, Opts(num_burn = burn_in, num_save = 1, update_sigma = FALSE,
                                 cache_trees = FALSE, num_print = .Machine$integer.max), warn = FALSE)
  next_stream(st)

  lambda0 <- 0.1; Rv <- rep(0, N); Wv <- exp(Rv); s1 <- 1; rho1 <- 0; H1 <- 0.25   # H1: rho_1 window
  if (seed != 20260816) {   # an additional chain starts from dispersed values
    lambda0 <- rgamma(1, 2, 10); Rv <- rnorm(N, 0, 0.5); Wv <- exp(Rv)
    s1 <- 1/rgamma(1, 2, 1); rho1 <- runif(1, -0.5, 0.5)
    if (spatial == 0) rho1 <- 0 }   # rho_1 is not a parameter of the independent model
  keep_lam <- keep_s1 <- keep_rho1 <- numeric(num_iter)
  acc_R <- rep(0, N); acc_rho1 <- 0
  d_cnt <- tabulate(cty[sd_df$delta == 1], nbins = N)
  sy    <- as.numeric(tapply(sd_df$time, factor(cty, levels = 1:N), sum)); sy[is.na(sy)] <- 0

  ## evaluation grid for S(t), DMS and the residuals
  tg  <- seq(t_hor/100, t_hor, length.out = 100); tg0 <- c(0, tg); nT0 <- length(tg0)
  dt  <- tg[1]
  cumtrap <- function(P, h) { m <- (P[, -1, drop = FALSE] + P[, -ncol(P), drop = FALSE])/2
                              t(apply(m, 1, cumsum)) * h }
  clip01 <- function(m) { m[m < 0] <- 0; m[m > 1] <- 1; m }
  ## int_0^a d(t) dt by the trapezoid rule on the grid dt, 2dt, ..., a, with d(0) = 0
  trap <- function(d) (sum(d) - d[length(d)]/2) * dt

  ## representative patient in the county of interest: the stratum's own medians
  in_c <- which(cty == county_oi)
  if (length(in_c) == 0) { lg(sprintf("WARNING: county %d has no subjects; using stratum medians", county_oi))
                           in_c <- seq_len(nt) }
  ## sample medians (lower median for the binary covariates and the grade)
  lowmed <- function(v) sort(v)[ceiling(length(v)/2)]
  rep_x <- function(rows) {
    Xs <- X[rows, , drop = FALSE]
    xr <- apply(Xs, 2, lowmed); xr["Age"] <- median(Xs[, "Age"])
    g  <- lowmed(1 + Xs[, "Grade2"] + 2*Xs[, "Grade3"])
    xr["Grade2"] <- as.numeric(g == 2); xr["Grade3"] <- as.numeric(g == 3)
    xr }
  x_rep <- rep_x(in_c)
  ## contrasts for DMS(a): flip one covariate at a time (TD, BD), and SMU best-vs-own
  x_TD <- x_rep; x_TD["TD"] <- 1;  x_TD0 <- x_rep; x_TD0["TD"] <- 0
  x_BD <- x_rep; x_BD["BD"] <- 1;  x_BD0 <- x_rep; x_BD0["BD"] <- 0
  M_best <- max(M_hat)
  ## prediction blocks: rows = scenarios, evaluated on tg0
  mk_block <- function(xv, Mv) cbind(rep(tg0, each = 1)/max_time, Mv,
                                     matrix(rep(xv, nT0), nrow = nT0, byrow = TRUE))
  SCEN <- list(rep_Mhat = list(x = x_rep,  M = M_hat[county_oi]),
               TD1      = list(x = x_TD,   M = M_hat[county_oi]),
               TD0      = list(x = x_TD0,  M = M_hat[county_oi]),
               BD1      = list(x = x_BD,   M = M_hat[county_oi]),
               BD0      = list(x = x_BD0,  M = M_hat[county_oi]),
               SMUbest  = list(x = x_rep,  M = M_best))
  scen_block <- do.call(rbind, lapply(SCEN, function(s) mk_block(s$x, s$M)))
  nS <- length(SCEN)

  ## additional prediction scenarios for the tables and figures (do not affect the MCMC)
  ## (a) Broward (county_oi) grade-specific contrasts, for Table 2
  xg <- function(xv, g) { xv["Grade2"] <- as.numeric(g == 2); xv["Grade3"] <- as.numeric(g == 3); xv }
  EXT <- list()
  for (g in 1:3) {
    xb <- xg(x_rep, g)
    EXT[[sprintf("g%d_rep", g)]] <- list(x = xb, M = M_hat[county_oi])
    EXT[[sprintf("g%d_SMU", g)]] <- list(x = xb, M = M_best)
    x1 <- xb; x1["TD"] <- 1; x0 <- xb; x0["TD"] <- 0
    EXT[[sprintf("g%d_TD1", g)]] <- list(x = x1, M = M_hat[county_oi])
    EXT[[sprintf("g%d_TD0", g)]] <- list(x = x0, M = M_hat[county_oi])
    x1 <- xb; x1["BD"] <- 1; x0 <- xb; x0["BD"] <- 0
    EXT[[sprintf("g%d_BD1", g)]] <- list(x = x1, M = M_hat[county_oi])
    EXT[[sprintf("g%d_BD0", g)]] <- list(x = x0, M = M_hat[county_oi])
  }
  ext_block <- do.call(rbind, lapply(EXT, function(s) mk_block(s$x, s$M)))
  nE <- length(EXT)
  ## (b) per-county representative-patient contrasts (every county with >= 1 subject)
  n_i <- tabulate(cty, nbins = N)
  cty_list <- which(n_i > 0)
  CS <- list(); cs_cty <- integer(0)
  for (cc in cty_list) {
    xc <- rep_x(which(cty == cc))
    x1 <- xc; x1["TD"] <- 1; x0 <- xc; x0["TD"] <- 0
    y1 <- xc; y1["BD"] <- 1; y0 <- xc; y0["BD"] <- 0
    CS <- c(CS, list(list(x = xc, M = M_hat[cc]), list(x = xc, M = M_best),
                     list(x = x1, M = M_hat[cc]), list(x = x0, M = M_hat[cc]),
                     list(x = y1, M = M_hat[cc]), list(x = y0, M = M_hat[cc])))
    cs_cty <- c(cs_cty, rep(cc, 6))
  }
  cs_block <- do.call(rbind, lapply(CS, function(s) mk_block(s$x, s$M)))
  nC <- length(cty_list)
  lg(sprintf("prediction scenarios: %d by grade, %d counties x 6", nE, nC))

  ## Step-3 quadrature design (also reused for the Cox-Snell residuals)
  qw <- (seq_len(NQ) - 0.5)/NQ
  q_rows  <- rep(1:nt, times = NQ)
  Xq_base <- cbind(as.numeric(outer(sd_df$time, qw))/max_time, 0, X[q_rows, , drop = FALSE])
  Xy_base <- cbind(sd_df$time/max_time, 0, X)

  n_keep_max <- sum(seq_len(num_iter) > burn_in & seq_len(num_iter) %% THIN == 0)
  S_scen_sum <- matrix(0, nS, length(tg)); S_scenM_sum <- matrix(0, nS, length(tg))
  vip_sum <- rep(0, ncol(Xd)); names(vip_sum) <- colnames(Xd)
  csnell_sum <- rep(0, nt)
  S_scen_draws <- array(NA_real_, c(nS, length(tg), n_keep_max))
  lw_keep <- numeric(n_keep_max); n_keep <- 0
  sumw <- 0; maxlw <- -Inf
  w_lam <- 0; w_W <- rep(0, N); w_rho1 <- 0; w_s1 <- 0; w_M <- rep(0, N)
  ext_draws <- array(NA_real_, c(nE, length(tg), n_keep_max))
  cty_dms <- array(NA_real_, c(nC, 6, n_keep_max),
                   dimnames = list(cty_list, c("SMU10","TD10","BD10","SMU5","TD5","BD5"), NULL))
  W_draws <- matrix(NA_real_, n_keep_max, N); M_draws <- matrix(NA_real_, n_keep_max, N)
  cty_surv <- array(NA_real_, c(nC, 2, n_keep_max), dimnames = list(cty_list, c("S5", "S10"), NULL))
  vip_draws <- matrix(NA_real_, n_keep_max, ncol(Xd)); colnames(vip_draws) <- colnames(Xd)
  R_trace <- matrix(NA_real_, num_iter, N)
  ia10 <- which(tg <= 10 + 1e-9); ia5 <- which(tg <= 5 + 1e-9)
  rescale_all <- function(f) {
    S_scen_sum <<- S_scen_sum*f; S_scenM_sum <<- S_scenM_sum*f; vip_sum <<- vip_sum*f
    csnell_sum <<- csnell_sum*f; w_lam <<- w_lam*f; w_W <<- w_W*f
    w_rho1 <<- w_rho1*f; w_s1 <<- w_s1*f; w_M <<- w_M*f; sumw <<- sumw*f }

  for (iter in 1:num_iter) {
    ## (a) data augmentation: rejected points of the thinned Poisson process
    rates <- lambda0 * Wv[cty] * sd_df$time
    q <- rpois(nt, pmax(rates, 1e-8))
    idx <- which(q > 0); Xaug <- NULL; m_cnt <- rep(0, N)
    if (length(idx) > 0) {
      rep_i <- rep(idx, q[idx]); G <- runif(length(rep_i), 0, sd_df$time[rep_i])
      XG <- cbind(G/max_time, M_hat[cty[rep_i]], X[rep_i, , drop = FALSE])
      keepG <- runif(nrow(XG)) > pnorm(forest$do_predict(XG))
      if (any(keepG)) { Xaug <- XG[keepG, , drop = FALSE]
                        m_cnt <- tabulate(cty[rep_i[keepG]], nbins = N) }
    }
    ## (b) probit back-fitting on events (1) and rejected points (0)
    ev_idx <- which(sd_df$delta == 1)
    Xall <- rbind(Xaug, Xd[ev_idx, , drop = FALSE])
    Yall <- c(rep(0, if (is.null(Xaug)) 0 else nrow(Xaug)), rep(1, length(ev_idx)))
    bcur <- forest$do_predict(Xall)
    Z <- ifelse(Yall == 1, rtruncnorm(length(Yall), a = 0, b = Inf, mean = bcur, sd = 1),
                           rtruncnorm(length(Yall), a = -Inf, b = 0, mean = bcur, sd = 1))
    st <- .Random.seed                 # do_gibbs uses R's RNG without saving its state:
    forest$do_gibbs(Xall, Z, Xall, 1)  # continue on the next independent stream
    next_stream(st)
    ## (c) closed-form updates
    lambda0 <- rgamma(1, 1 + sum(sd_df$delta) + sum(m_cnt), 1 + sum(Wv[cty]*sd_df$time))
    qf <- if (spatial == 1) as.numeric(t(Rv) %*% (D - rho1*A) %*% Rv) else sum(Rv^2)
    s1 <- 1/rgamma(1, 1 + N/2, 1 + qf/2)
    ## (d) rho_1 (only defined for the spatial model) and R = log W
    if (spatial == 1) {
      rp <- rho_prop(rho1, H1)
      if (log(runif(1)) < car_ld(Rv, rp, s1) - car_ld(Rv, rho1, s1) +
                          log(rho_win(rho1, H1)) - log(rho_win(rp, H1))) {
        rho1 <- rp; acc_rho1 <- acc_rho1 + 1 }
    }
    for (c2 in 1:N) {
      Rp <- Rv; Rp[c2] <- rnorm(1, Rv[c2], 0.2)
      ## every rejected point also carries a factor lambda_0 W_i
      llp <- (d_cnt[c2] + m_cnt[c2])*Rp[c2] - lambda0*exp(Rp[c2])*sy[c2]
      llc <- (d_cnt[c2] + m_cnt[c2])*Rv[c2] - lambda0*exp(Rv[c2])*sy[c2]
      la <- (llp + frail_ld(Rp, rho1, s1)) - (llc + frail_ld(Rv, rho1, s1))
      if (is.finite(la) && log(runif(1)) < la) { Rv <- Rp; acc_R[c2] <- acc_R[c2] + 1 }
    }
    Wv <- exp(Rv)          # no deterministic recentring: the CAR here is proper
    keep_lam[iter] <- lambda0; keep_s1[iter] <- s1; keep_rho1[iter] <- rho1
    R_trace[iter, ] <- Rv

    ## ==== 3. STEP 3 -- importance weight for this draw, and the weighted summaries
    if (iter > burn_in && iter %% THIN == 0) {
      lw <- 0; Mdraw <- M_hat
      if (use_step3 == 1) {
        Mdraw <- M_post[sample.int(nrow(M_post), 1), ]
        Xq_m <- Xq_base; Xq_m[, 2] <- Mdraw[cty[q_rows]]
        Xq_h <- Xq_base; Xq_h[, 2] <- M_hat[cty[q_rows]]
        Xy_m <- Xy_base; Xy_m[, 2] <- Mdraw[cty]
        Xy_h <- Xy_base; Xy_h[, 2] <- M_hat[cty]
        Im  <- rowMeans(matrix(pnorm(forest$do_predict(Xq_m)), nt, NQ)) * sd_df$time
        Ih  <- rowMeans(matrix(pnorm(forest$do_predict(Xq_h)), nt, NQ)) * sd_df$time
        lPm <- pnorm(forest$do_predict(Xy_m), log.p = TRUE)
        lPh <- pnorm(forest$do_predict(Xy_h), log.p = TRUE)
        lw  <- sum(sd_df$delta*(lPm - lPh)) - lambda0*sum(Wv[cty]*(Im - Ih))
      } else {
        Ih <- rowMeans(matrix(pnorm(forest$do_predict({z <- Xq_base; z[,2] <- M_hat[cty[q_rows]]; z})),
                              nt, NQ)) * sd_df$time
      }
      if (lw > maxlw) { if (is.finite(maxlw)) rescale_all(exp(maxlw - lw)); maxlw <- lw; w <- 1
      } else w <- exp(lw - maxlw)
      sumw <- sumw + w
      w_lam <- w_lam + w*lambda0; w_W <- w_W + w*Wv
      w_rho1 <- w_rho1 + w*rho1;  w_s1 <- w_s1 + w*s1; w_M <- w_M + w*Mdraw
      ## variable inclusion proportions (weighted, after burn-in)
      cnt <- as.numeric(forest$get_counts())
      vip_sum <- vip_sum + w * (if (sum(cnt) > 0) cnt/sum(cnt) else cnt)
      ## Cox-Snell residual r_ij = Lambda_hat(y_ij) = lambda_0 W_i Int_0^y Phi(b)
      csnell_sum <- csnell_sum + w * (lambda0 * Wv[cty] * (if (use_step3 == 1) Im else Ih))
      ## survival curves for the scenario grid, at M_hat and at the drawn M
      phS <- matrix(pnorm(forest$do_predict(scen_block)), nT0, nS)
      S_s <- clip01(exp(-lambda0 * Wv[county_oi] * cumtrap(t(phS), dt)))
      S_scen_sum <- S_scen_sum + w * S_s
      n_keep <- n_keep + 1; S_scen_draws[, , n_keep] <- S_s; lw_keep[n_keep] <- lw
      phE <- matrix(pnorm(forest$do_predict(ext_block)), nT0, nE)
      ext_draws[, , n_keep] <- clip01(exp(-lambda0 * Wv[county_oi] * cumtrap(t(phE), dt)))
      phC <- matrix(pnorm(forest$do_predict(cs_block)), nT0, length(CS))
      SC  <- clip01(exp(-(cumtrap(t(phC), dt) * (lambda0 * Wv[cs_cty]))))
      for (k in seq_len(nC)) {
        r0 <- (k - 1) * 6
        cty_dms[k, "SMU10", n_keep] <- trap(SC[r0+2, ia10] - SC[r0+1, ia10])
        cty_dms[k, "TD10", n_keep] <- trap(SC[r0+3, ia10] - SC[r0+4, ia10])
        cty_dms[k, "BD10", n_keep] <- trap(SC[r0+5, ia10] - SC[r0+6, ia10])
        cty_dms[k, "SMU5", n_keep] <- trap(SC[r0+2, ia5] - SC[r0+1, ia5])
        cty_dms[k, "TD5", n_keep] <- trap(SC[r0+3, ia5] - SC[r0+4, ia5])
        cty_dms[k, "BD5", n_keep] <- trap(SC[r0+5, ia5] - SC[r0+6, ia5])
      }
      cty_surv[, "S5", n_keep]  <- SC[(seq_len(nC) - 1) * 6 + 1, max(ia5)]
      cty_surv[, "S10", n_keep] <- SC[(seq_len(nC) - 1) * 6 + 1, max(ia10)]
      W_draws[n_keep, ] <- Wv; M_draws[n_keep, ] <- Mdraw
      vip_draws[n_keep, ] <- cnt
      if (use_step3 == 1) {
        ## county scenarios take the drawn M_i; the best-county scenario stays fixed
        sb <- scen_block; sb[, 2] <- c(rep(Mdraw[county_oi], 5*nT0), rep(M_best, nT0))
        phM <- matrix(pnorm(forest$do_predict(sb)), nT0, nS)
        S_scenM_sum <- S_scenM_sum + w * clip01(exp(-lambda0*Wv[county_oi]*cumtrap(t(phM), dt)))
      }
    }
    if (iter %% 1000 == 0)
      lg(sprintf("iter %6d  lam0=%.3f s1=%.2f rho1=%.2f accR=%.2f (%.1f min)",
                 iter, lambda0, s1, rho1, mean(acc_R)/iter,
                 as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }

  ## ==== 4. Posterior summaries -- every one of them Step-3 weighted
  lwk <- lw_keep[1:n_keep]; Wt <- exp(lwk - max(lwk))
  ess_w <- sum(Wt)^2/sum(Wt^2)
  w_top <- max(Wt)/sum(Wt); p_zero <- mean(lwk == 0)
  if (use_step3 == 1) {
    lg(sprintf(paste("Step 3: %d retained draws, ESS = %.1f (%.0f%%), largest normalised weight %.3f,",
                     "share of log-weights = 0: %.2f"), n_keep, ess_w, 100*ess_w/n_keep, w_top, p_zero))
  } else lg(sprintf("Step 3 OFF (2-step method, unit weights): %d retained draws", n_keep))
  if (ess_w < 0.1*n_keep)
    lg("WARNING: importance weights are close to degenerate; treat interval estimates with caution")

  S_rep   <- S_scen_sum/sumw                     # posterior mean survival, per scenario
  S_repM  <- if (use_step3 == 1) S_scenM_sum/sumw else NULL
  vip_raw <- vip_sum/sumw
  ## VIP of t and the 6 covariates (grade = its two indicators combined)
  vip6 <- c(Age = unname(vip_raw["Age"]), HR = unname(vip_raw["HR"]),
            Grade = unname(vip_raw["Grade2"] + vip_raw["Grade3"]), SMU = unname(vip_raw["M"]),
            TD = unname(vip_raw["TD"]), BD = unname(vip_raw["BD"]))
  csnell <- csnell_sum/sumw
  devres <- dev_res(sd_df$delta, csnell)

  ## DMS(a) = Int_0^a [S(t | x*, M*) - S(t | x, M)] dt, per draw, then weighted
  dms_of <- function(k_star, k_ref, a) {
    ia <- which(tg <= a)
    apply(S_scen_draws[, , 1:n_keep, drop = FALSE], 3, function(S)
      trap(S[k_star, ia] - S[k_ref, ia])) }
  sn <- names(SCEN)
  dms_tab <- do.call(rbind, lapply(list(c("TD1","TD0"), c("BD1","BD0"), c("SMUbest","rep_Mhat")),
    function(pr) {
      d <- dms_of(match(pr[1], sn), match(pr[2], sn), t_hor)
      h <- whpd(d, Wt)
      data.frame(contrast = paste(pr[1], "vs", pr[2]),
                 DMS = sum(Wt*d)/sum(Wt), HPD_low = h[1], HPD_high = h[2]) }))
  lg(sprintf("DMS(%g): %s", t_hor,
             paste(sprintf("%s = %.3f (%.3f, %.3f)", dms_tab$contrast, dms_tab$DMS,
                           dms_tab$HPD_low, dms_tab$HPD_high), collapse = " | ")))
  ## pointwise 95% credible band for the representative patient's survival curve
  qq <- apply(S_scen_draws[1, , 1:n_keep, drop = FALSE], 2, function(v) whpd(as.numeric(v), Wt))
  S_rep_lo <- qq[1, ]; S_rep_hi <- qq[2, ]

  saveRDS(list(
    race = race, stage = stage, n = nt, events = sum(sd_df$delta),
    settings = c(num_iter = num_iter, burn_in = burn_in, num_tree = num_tree,
                 THIN = THIN, NQ = NQ, n1 = n1, b1 = b1, spatial = spatial,
                 use_step3 = use_step3, county = county_oi, t_hor = t_hor, seed = seed),
    n_steps = if (use_step3 == 1) 3L else 2L,
    n_keep = n_keep, ess_weights = ess_w, w_top = w_top, p_logw_zero = p_zero, logw = lwk,
    M_hat = M_hat, M_post_mean_step3 = w_M/sumw, var_logitM = var_logitM, n_0i = n_0i,
    W_hat = w_W/sumw, lam0_post = w_lam/sumw, s1_post = w_s1/sumw, rho1_post = w_rho1/sumw,
    vip_raw = vip_raw, vip6 = vip6,
    times = tg, S_rep = S_rep, S_rep_Mdrawn = S_repM,
    S_rep_lo = S_rep_lo, S_rep_hi = S_rep_hi, scenarios = sn,
    dms = dms_tab, cox_snell = csnell, deviance_res = devres,
    acc_R = mean(acc_R)/num_iter, acc_rho1 = acc_rho1/num_iter,
    minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
    file.path(outdir, sprintf("fcr_%s_%s.rds", race, stage)))

  lg(sprintf("DONE %s-%s (%d-step, %s frailty) | ESS=%.1f | VIP top: %s | %.1f min",
             race, stage, if (use_step3 == 1) 3L else 2L,
             if (spatial == 1) "spatial" else "independent", ess_w,
             paste(names(sort(vip6, decreasing = TRUE))[1:2], collapse = ","),
             as.numeric(difftime(Sys.time(), t_start, units = "mins"))))


  ## ==== 5. Detailed outputs for report (also DMS(5) and unweighted summaries)
  summ <- function(v, w) { h <- whpd(v, w)
    c(mean = sum(w*v)/sum(w), hpd_lo = h[1], hpd_hi = h[2], q025 = wq(v, w, .025), q975 = wq(v, w, .975)) }
  U <- rep(1, n_keep)
  en <- names(EXT); exS <- ext_draws[, , 1:n_keep, drop = FALSE]
  dms_e <- function(ia, k1, k0) apply(exS, 3, function(S) trap(S[k1, ia] - S[k0, ia]))
  tab2 <- do.call(rbind, lapply(1:3, function(g) do.call(rbind, lapply(list(
      c("SMU", sprintf("g%d_SMU", g), sprintf("g%d_rep", g)),
      c("BD",  sprintf("g%d_BD1", g), sprintf("g%d_BD0", g)),
      c("TD",  sprintf("g%d_TD1", g), sprintf("g%d_TD0", g))), function(pr) {
    k1 <- match(pr[2], en); k0 <- match(pr[3], en)
    d10 <- dms_e(ia10, k1, k0); d5 <- dms_e(ia5, k1, k0)
    s10 <- summ(d10, Wt); s5 <- summ(d5, Wt); u10 <- summ(d10, U); u5 <- summ(d5, U)
    data.frame(stratum = paste(race, stage), grade = g, contrast = pr[1],
               DMS10 = s10["mean"], HPD10_lo = s10["hpd_lo"], HPD10_hi = s10["hpd_hi"],
               DMS5 = s5["mean"], HPD5_lo = s5["hpd_lo"], HPD5_hi = s5["hpd_hi"],
               DMS10_unweighted = u10["mean"], HPD10_lo_unw = u10["hpd_lo"], HPD10_hi_unw = u10["hpd_hi"],
               DMS5_unweighted = u5["mean"], row.names = NULL) }))))
  ## per-draw values of the same contrasts
  cn_g <- expand.grid(contrast = c("SMU", "BD", "TD"), grade = 1:3, stringsAsFactors = FALSE)
  k_of <- function(cn, g, star) match(sprintf(if (cn == "SMU") (if (star) "g%d_SMU" else "g%d_rep")
                                            else paste0("g%d_", cn, if (star) "1" else "0"), g), en)
  dms_draws <- function(ia) sapply(seq_len(nrow(cn_g)), function(r)
    dms_e(ia, k_of(cn_g$contrast[r], cn_g$grade[r], TRUE), k_of(cn_g$contrast[r], cn_g$grade[r], FALSE)))
  tab2_draws <- list(DMS10 = dms_draws(ia10), DMS5 = dms_draws(ia5))
  colnames(tab2_draws$DMS10) <- colnames(tab2_draws$DMS5) <- paste0(cn_g$contrast, "_g", cn_g$grade)
  cd <- cty_dms[, , 1:n_keep, drop = FALSE]
  cty_tab <- do.call(rbind, lapply(seq_len(nC), function(k) {
    row <- data.frame(county = cty_list[k], n_i = n_i[cty_list[k]], M_hat = M_hat[cty_list[k]],
                      n_0i = n_0i[cty_list[k]])
    for (nm in dimnames(cd)[[2]]) { v <- cd[k, nm, ]; s <- summ(v, Wt); u <- summ(v, U)
      row[[nm]] <- s["mean"]; row[[paste0(nm, "_hpd_lo")]] <- s["hpd_lo"]; row[[paste0(nm, "_hpd_hi")]] <- s["hpd_hi"]
      row[[paste0(nm, "_q025")]] <- s["q025"]; row[[paste0(nm, "_q975")]] <- s["q975"]
      row[[paste0(nm, "_unw")]] <- u["mean"] }
    row }))
  Wd <- W_draws[1:n_keep, , drop = FALSE]
  W_tab <- data.frame(county = 1:N, n_i = n_i,
    W_mean = apply(Wd, 2, function(v) sum(Wt*v)/sum(Wt)),
    W_q025 = apply(Wd, 2, function(v) wq(v, Wt, .025)), W_q975 = apply(Wd, 2, function(v) wq(v, Wt, .975)),
    W_mean_unw = colMeans(Wd), W_q025_unw = apply(Wd, 2, function(v) wq(v, U, .025)),
    W_q975_unw = apply(Wd, 2, function(v) wq(v, U, .975)))
  vd <- vip_draws[1:n_keep, , drop = FALSE]
  saveRDS(list(race = race, stage = stage, n = nt, events = sum(sd_df$delta),
    settings = c(num_iter = num_iter, burn_in = burn_in, num_tree = num_tree, THIN = THIN, NQ = NQ,
                 n1 = n1, b1 = b1, spatial = spatial, use_step3 = use_step3, county = county_oi, t_hor = t_hor,
                 seed = seed),
    ess_weights = ess_w, w_top = w_top, p_logw_zero = p_zero, n_keep = n_keep, logw = lwk,
    tab2 = tab2, cty_tab = cty_tab, W_tab = W_tab,
    tab2_draws = tab2_draws, cty_dms_draws = cd, cty_surv_draws = cty_surv[, , 1:n_keep, drop = FALSE],
    S_rep_draws = S_scen_draws[1, , 1:n_keep], cty_list = cty_list,
    vip6_weighted = vip_of(vd, Wt), vip6_unweighted = vip_of(vd, U), vip_draws = vd,
    S_rep = S_rep, S_rep_lo = S_rep_lo, S_rep_hi = S_rep_hi, times = tg, scenarios = sn, dms = dms_tab,
    S_rep_unweighted = apply(S_scen_draws[, , 1:n_keep, drop = FALSE], c(1, 2), mean),
    cox_snell = csnell, deviance_res = devres, delta = sd_df$delta, county = cty, Xraw = Xraw,
    M_hat = M_hat, var_logitM = var_logitM, n_0i = n_0i, M_post_q = apply(M_post, 2, quantile, c(.025, .5, .975)),
    lam_trace = keep_lam, s1_trace = keep_s1, rho1_trace = keep_rho1, R_trace = R_trace,
    rho0_trace = rho0_tr, s0_trace = s0_tr,
    W_draws = Wd, M_draws = M_draws[1:n_keep, , drop = FALSE], x_rep = x_rep, M_best = M_best,
    acc_R = mean(acc_R)/num_iter, acc_rho1 = acc_rho1/num_iter,
    acc_logitM = mean(acc1)/n1, acc_rho0 = acc_r0/n1,
    minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
    file.path(outdir, sprintf("fcr_%s_%s_details.rds", race, stage)))
  lg("detailed outputs saved")
}

################################################################################
# PART B -- report: pooled posterior summaries, tables, figures and numbers
# Chains share Step 1, so the pooled weights are exp(logw - max over all draws).
################################################################################
make_report <- function(runs_dir, figdir) {
  suppressMessages({ library(survival); library(readxl); library(maps) })
  dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
  STRATA <- c("AA_L", "WA_L", "AA_R", "WA_R", "AA_D", "WA_D")
  STAGE <- c(L = "Local", R = "Regional", D = "Distant")
  N <- 67
  COUNTY <- c("Alachua","Baker","Bay","Bradford","Brevard","Broward","Calhoun","Charlotte","Citrus","Clay",
              "Collier","Columbia","Miami-Dade","DeSoto","Dixie","Duval","Escambia","Flagler","Franklin","Gadsden",
              "Gilchrist","Glades","Gulf","Hamilton","Hardee","Hendry","Hernando","Highlands","Hillsborough","Holmes",
              "Indian River","Jackson","Jefferson","Lafayette","Lake","Lee","Leon","Levy","Liberty","Madison",
              "Manatee","Marion","Martin","Monroe","Nassau","Okaloosa","Okeechobee","Orange","Osceola","Palm Beach",
              "Pasco","Pinellas","Polk","Putnam","St. Johns","St. Lucie","Santa Rosa","Sarasota","Seminole","Sumter",
              "Suwannee","Taylor","Union","Volusia","Wakulla","Walton","Washington")   # county codes 1-67 (FIPS order)
  BROWARD <- 6                                                          # the county of interest of the runs
  MARK <- c(Baker = 2, Escambia = 17, Broward = 6, Citrus = 9)          # counties marked in Figures 5 and S11
  MARK_PCH <- c(Baker = 21, Escambia = 22, Broward = 8, Citrus = 24)
  LONG_DELAY <- 3                                                       # data code of a long TD / BD
  ## quartz embeds the fonts
  fig <- function(f, w, h) {
    if (capabilities("aqua")) quartz(type = "pdf", file = file.path(figdir, f), width = w, height = h)
    else grDevices::pdf(file.path(figdir, f), width = w, height = h) }
  wm <- function(v, w) sum(w*v)/sum(w)
  f3 <- function(v) sub("^-(0\\.0+)$", "\\1", sprintf("%.3f", v))     # 3 decimals, never "-0.000"
  est <- function(v, w) { h <- whpd(v, w); sprintf("%s (%s, %s)", f3(wm(v, w)), f3(h[1]), f3(h[2])) }
  rep_grade <- function(P) if (P$x_rep[["Grade3"]] == 1) 3 else if (P$x_rep[["Grade2"]] == 1) 2 else 1

  ## ---- read and pool the chains
  load_stratum <- function(S) {
    f <- function(d) file.path(runs_dir, d, sprintf("fcr_%s_details.rds", S))
    dirs <- list.files(runs_dir, pattern = sprintf("^%s(_chain[0-9]+)?$", S))
    dirs <- dirs[file.exists(f(dirs))]
    if (!length(dirs)) return(NULL)
    k <- suppressWarnings(as.integer(sub("^.*_chain", "", dirs)))
    k[is.na(k)] <- 1L
    dirs <- dirs[order(k)]                                              # S, S_chain2, S_chain3, ...
    X <- lapply(f(dirs), readRDS)
    names(X) <- dirs
    ## the chains may differ only in their length and seed
    same <- c("n1", "b1", "num_tree", "NQ", "spatial", "use_step3", "county", "t_hor")
    for (x in X)
      stopifnot("report needs runs with Step 3 (use_step3 = 1)" = x$settings[["use_step3"]] == 1,
                "report needs runs for the county of interest 6 (Broward)" = x$settings[["county"]] == BROWARD,
                "the chains of a stratum must share Step 1 and the model" =
                  identical(x$settings[same], X[[1]]$settings[same]))
    seeds <- sapply(X, function(x) x$settings[["seed"]])
    stopifnot("each chain of a stratum needs its own seed" = !anyDuplicated(seeds))
    for (x in X[-1])
      stopifnot(identical(x$delta, X[[1]]$delta), isTRUE(all.equal(x$M_hat, X[[1]]$M_hat)),
                identical(x$cty_list, X[[1]]$cty_list), isTRUE(all.equal(x$x_rep, X[[1]]$x_rep)),
                isTRUE(all.equal(x$times, X[[1]]$times)))
    lw <- unlist(lapply(X, `[[`, "logw"))
    w <- exp(lw - max(lw))
    ch <- rep(seq_along(X), sapply(X, function(x) length(x$logw)))
    Tc <- sapply(split(w, ch), sum)                                     # each chain's share of the pooled weight
    ## Cox-Snell residuals: each chain's weighted mean, combined with the chain's total weight
    cs <- Reduce(`+`, Map(function(x, t) t * x$cox_snell, X, Tc))/sum(Tc)
    cd <- lapply(X, `[[`, "cty_dms_draws")
    list(X = X, w = w, n = X[[1]]$n, x_rep = X[[1]]$x_rep, M_hat = X[[1]]$M_hat, times = X[[1]]$times,
         delta = X[[1]]$delta, county = X[[1]]$county, cty_list = X[[1]]$cty_list,
         vip_draws = do.call(rbind, lapply(X, `[[`, "vip_draws")),
         DMS10 = do.call(rbind, lapply(X, function(x) x$tab2_draws$DMS10)),
         S_rep_draws = do.call(cbind, lapply(X, `[[`, "S_rep_draws")),
         cty_dms = array(unlist(cd), c(dim(cd[[1]])[1:2], length(w)),
                         dimnames = c(dimnames(cd[[1]])[1:2], list(NULL))),
         W_draws = do.call(rbind, lapply(X, `[[`, "W_draws")),
         cox_snell = cs, deviance_res = dev_res(X[[1]]$delta, cs))
  }
  R <- setNames(lapply(STRATA, load_stratum), STRATA)
  have <- STRATA[!sapply(R, is.null)]
  if (!length(have)) stop("no runs found in ", runs_dir)
  for (s in have) cat(s, "pools", paste(names(R[[s]]$X), collapse = ", "), "\n")
  miss <- setdiff(STRATA, have)
  if (length(miss))
    message("Warning: no runs for ", paste(miss, collapse = ", "), " in ", runs_dir, ": skipped, so the ",
            "tables, figures and numbers that need these strata are incomplete or not written")

  ## DMS(10) of the representative patient of every county with patients
  county_table <- function(P) {
    n_i <- tabulate(P$county, nbins = N)
    do.call(rbind, lapply(seq_along(P$cty_list), function(k) {
      cc <- P$cty_list[k]
      r <- data.frame(county = cc, name = COUNTY[cc], n_i = n_i[cc])
      for (v in c("SMU10", "TD10", "BD10")) {
        d <- P$cty_dms[k, v, ]
        h <- whpd(d, P$w)
        r[[v]] <- wm(d, P$w)
        r[[paste0(v, "_lo")]] <- h[1]
        r[[paste0(v, "_hi")]] <- h[2]
      }
      r }))
  }
  CT <- lapply(R[have], county_table)

  ## ---- data (as in fit) and the competing PH-frailty model
  d0 <- read_fcr()
  d0$y <- pmax(d0$days/365.25, 1/365.25)
  d0$race <- factor(d0$Race, levels = c(2, 1), labels = c("AA", "WA"))    # other codes are NA, as in fit
  d0$stage <- factor(STAGE[d0$Stage], STAGE)
  stratum_data <- function(S, P) {
    race_code <- if (substr(S, 1, 2) == "AA") 2 else 1
    ds <- d0[d0$Race == race_code & d0$Stage == match(substr(S, 4, 4), names(STAGE)), ]
    ds <- ds[!is.na(ds$county) & ds$county >= 1 & ds$county <= N & !is.na(ds$days) & !is.na(ds$death), ]
    stopifnot("the data do not match the runs" = nrow(ds) == length(P$delta) &&
                all((ds$death == 1) == (P$delta == 1)) && all(ds$county == P$county))
    ds$Mhat <- P$M_hat[ds$county]
    ds }
  ## PH-frailty: Cox model with M_hat and iid Gaussian county frailties; cox.zph is the
  ## Grambsch-Therneau test of proportional hazards
  PH <- lapply(setNames(have, have), function(S) {
    ds <- stratum_data(S, R[[S]])
    ph <- coxph(Surv(y, death == 1) ~ Age + I(HR == 1) + I(Tg == 2) + I(Tg == 3) + I(TX == 1) + I(BX == 1) +
                  Mhat + frailty(county, distribution = "gaussian"),
                data = ds, control = coxph.control(iter.max = 1000))
    list(ds = ds, cf = suppressWarnings(summary(ph))$coefficients,    # (its unused column se2 can be NaN)
         zph = suppressWarnings(cox.zph(ph))$table,
         dev = residuals(ph, type = "deviance"),
         cs = (ds$death == 1) - residuals(ph, type = "martingale")) })
  ph_hr <- function(S, term) {
    cf <- PH[[S]]$cf
    sprintf("%.2f (%.2f, %.2f), p = %.2g", exp(cf[term, 1]), exp(cf[term, 1] - 1.96*cf[term, 2]),
            exp(cf[term, 1] + 1.96*cf[term, 2]), cf[term, "p"]) }

  ## ==== Main text: Table 1, Figures 3-5, Table 2
  ## ---- Table 1: VIP of the time t and the 6 covariates
  VIP <- sapply(STRATA, function(s) if (s %in% have) vip_of(R[[s]]$vip_draws, R[[s]]$w) else rep(NA, 7))
  rownames(VIP) <- c("Age", "HR", "Grade", "SMU", "TD", "BD", "time")
  t1 <- c(time = "Time ($t$)", Age = "Age", HR = "HR Status", Grade = "Tumor Grade", SMU = "SMU",
          TD = "Treatment Delay", BD = "Biopsy Delay")
  write.csv(data.frame(variable = names(t1), round(VIP[names(t1), ], 3), check.names = FALSE),
            file.path(figdir, "table1_vip.csv"), row.names = FALSE)
  cell <- function(v) ifelse(is.na(v), "--", f3(v))
  head1 <- paste(sprintf("(%s, %s)", substr(STRATA, 1, 2), substr(STRATA, 4, 4)), collapse = " & ")
  row1 <- function(k) paste0(t1[k], " & \\vline & ", paste(cell(VIP[k, ]), collapse = " & "), "\\\\")
  writeLines(c("% Table 1: posterior VIP (Step-3 weighted, chains pooled)",
               "\\begin{tabular}{cccccccc}", "\\hline",
               paste0("Strata: (Race, Stage) & & ", head1, "\\\\"), "\\hline",
               sapply(names(t1), row1),
               "\\hline", "\\end{tabular}"), file.path(figdir, "table1_vip.tex"))

  ## ---- Figure 3: county proportions of long BD, No SMU (1 - M_hat) and long TD
  Mhat_race <- function(r) {
    s <- intersect(paste0(r, c("_L", "_R", "_D")), have)
    if (length(s)) R[[s[1]]]$M_hat else numeric(0) }
  vals <- list()
  for (r in c("AA", "WA")) {
    s <- d0[d0$Race == (if (r == "AA") 2 else 1), ]
    vals[[r]] <- list(`Biopsy Delay` = tapply(s$BX == LONG_DELAY, s$county, mean),
                      `No SMU` = 1 - Mhat_race(r),
                      `Treatment Delay` = tapply(s$TX == LONG_DELAY, s$county, mean))
  }
  fig("Figure3_boxplots_3step.pdf", 6.2, 3.5)
  par(mfrow = c(1, 3), mar = c(3, 3.4, 2.2, 0.6), mgp = c(2.2, 0.7, 0), las = 1)
  par(cex = 0.85)
  for (v in c("Biopsy Delay", "No SMU", "Treatment Delay")) {
    boxplot(list(AA = vals$AA[[v]], WA = vals$WA[[v]]), ylim = c(0, 1), main = v, cex.main = 1.1,
            col = "white", pch = 16, cex = 0.6)
    if (v == "Biopsy Delay") title(ylab = "Proportion")
  }
  invisible(dev.off())

  ## ---- Figure 4: posterior mean survival of the representative patients in Broward
  S_mean <- lapply(R[have], function(P) colSums(P$w * t(P$S_rep_draws))/sum(P$w))
  fig("Figure4_survival_3step.pdf", 7, 3.1)
  par(mfrow = c(1, 3), mar = c(3.6, 4, 2, 0.6), mgp = c(2.3, 0.7, 0), las = 1)
  par(cex = 0.85)
  for (st in names(STAGE)) {
    plot(NA, xlim = c(0, 10), ylim = c(0, 1), xlab = "Years",
         ylab = if (st == "L") expression(hat(S)(t)) else "", main = STAGE[st])
    for (r in c("AA", "WA")) {
      S <- paste0(r, "_", st)
      if (!is.null(S_mean[[S]]))
        lines(c(0, R[[S]]$times), c(1, S_mean[[S]]), lty = if (r == "AA") 2 else 1, lwd = 1.6)
    }
    if (st == "L") legend("topright", c("AA", "WA"), lty = c(2, 1), lwd = 1.6, bty = "n")
  }
  invisible(dev.off())

  ## ---- Table 2: DMS(10) of the representative AA patient in Broward, by grade and stage
  t2 <- expand.grid(grade = 1:3, contrast = c("SMU", "BD", "TD"), stringsAsFactors = FALSE)[, 2:1]
  t2_row <- function(s, i) {
    if (!s %in% have) return(data.frame(stratum = s, t2[i, ], DMS10 = NA, HPD_lo = NA, HPD_hi = NA))
    v <- R[[s]]$DMS10[, sprintf("%s_g%d", t2$contrast[i], t2$grade[i])]
    h <- whpd(v, R[[s]]$w)
    data.frame(stratum = s, t2[i, ], DMS10 = wm(v, R[[s]]$w), HPD_lo = h[1], HPD_hi = h[2]) }
  T2 <- do.call(rbind, lapply(paste0("AA_", names(STAGE)), function(s)
    do.call(rbind, lapply(seq_len(nrow(t2)), function(i) t2_row(s, i)))))
  T2$stage <- STAGE[substr(T2$stratum, 4, 4)]
  write.csv(data.frame(T2[, c("stage", "contrast", "grade")], round(T2[, c("DMS10", "HPD_lo", "HPD_hi")], 3)),
            file.path(figdir, "table2_dms10_broward.csv"), row.names = FALSE)
  lab2 <- list(SMU = c("SMU:", "$M_i^*=M^{\\max},$", "$M_i=\\widehat{M}_i$"),
               BD  = c("Biopsy Delay:", "Short vs Long;", "$M_i^*=M_i=\\widehat{M}_i$"),
               TD  = c("Treatment Delay:", "Short vs Long;", "$M_i^*=M_i=\\widehat{M}_i$"))
  row2 <- function(cn, g) {
    x <- T2[T2$contrast == cn & T2$grade == g, ]
    cells <- ifelse(is.na(x$DMS10), "--", sprintf("$%s$ $(%s, %s)$", f3(x$DMS10), f3(x$HPD_lo), f3(x$HPD_hi)))
    paste0(lab2[[cn]][g], " & ", g, " & ", paste(cells, collapse = " & "), " \\\\") }
  writeLines(c(paste("% Table 2: DMS(10) in years (95% HPD), representative AA patient in Broward",
                     "(Step-3 weighted, chains pooled)"),
               "\\begin{tabular}{l|c|c|c|c}", "\\toprule",
               "Covariate with & Tumor & \\multicolumn{3}{c}{\\underline{Cancer Stage}} \\\\",
               "different values & Grade & Local & Regional & Distant \\\\",
               unlist(lapply(c("SMU", "BD", "TD"), function(cn) c("\\midrule", sapply(1:3, function(g) row2(cn, g))))),
               "\\bottomrule", "\\end{tabular}"), file.path(figdir, "table2_dms10_broward.tex"))

  ## ---- Figures 5 and S11: county DMS(10) for TD against n_i
  plot_td_by_n <- function(file, ct, mar_r, ylim, ylab) {
    fig(file, 6.5, 4.1)
    par(mar = c(4.2, 4.4, 0.8, mar_r), las = 1)
    plot(ct$n_i, ct$TD10, type = "n", log = "x", ylim = ylim,
         xlab = expression("County size " * n[i] * " (log scale)"), ylab = ylab)
    abline(h = 0, lty = 3)
    segments(ct$n_i, ct$TD10_lo, ct$n_i, ct$TD10_hi, col = "grey55", lwd = 0.8)
    o <- !(ct$county %in% MARK)
    points(ct$n_i[o], ct$TD10[o], pch = 16, cex = 0.55)
    for (nm in names(MARK)) {                                           # the marked counties in black, on top
      j <- ct$county == MARK[nm]
      if (any(j)) {
        segments(ct$n_i[j], ct$TD10_lo[j], ct$n_i[j], ct$TD10_hi[j], col = "black", lwd = 1.3)
        points(ct$n_i[j], ct$TD10[j], pch = MARK_PCH[nm], bg = "white", cex = 1.2, lwd = 1.2)
      }
    }
    invisible(dev.off()) }
  if ("AA_L" %in% have) {
    ct <- CT$AA_L
    plot_td_by_n("Figure5_dms10_td_3step.pdf", ct, 1.8, range(0, ct$TD10_lo, ct$TD10_hi), "DMS(10) (years)")
  }

  ## ==== Supplement Section 2: Figures S5-S17
  ## ---- Figure S5: Kaplan-Meier curves by stage and race; log-rank tests
  fig("FigureS5_km.pdf", 7, 3.1)
  par(mfrow = c(1, 3), mar = c(3.6, 4, 2, 0.6), mgp = c(2.3, 0.7, 0), las = 1)
  par(cex = 0.85)
  LR <- list()
  for (s in levels(d0$stage)) {
    ds <- d0[d0$stage == s, ]
    km <- survfit(Surv(y, death == 1) ~ race, data = ds)
    plot(km, lty = c(2, 1), lwd = 1.4, xlim = c(0, 10), ylim = c(0, 1), xlab = "Years",
         ylab = if (s == "Local") expression(hat(S)[KM](t)) else "", main = s)
    if (s == "Local") legend("bottomleft", c("AA", "WA"), lty = c(2, 1), lwd = 1.4, bty = "n")
    LR[[s]] <- survdiff(Surv(y, death == 1) ~ race, data = ds)$chisq
  }
  invisible(dev.off())
  LR_strat <- survdiff(Surv(y, death == 1) ~ race + strata(stage), data = d0)$chisq

  ## ---- Figures S6-S8: residuals, spatial SBART vs PH-frailty
  RES <- list()
  for (S in intersect(c("AA_R", "WA_R"), have))
    RES[[S]] <- list(ds = PH[[S]]$ds, cs_sb = R[[S]]$cox_snell, dv_sb = R[[S]]$deviance_res,
                     cs_ph = PH[[S]]$cs, dv_ph = PH[[S]]$dev)
  if (length(RES)) {
    fig("FigureS6_coxsnell.pdf", 7, 3.6)
    par(mfrow = c(1, 2), mar = c(3.6, 4, 2, 0.6), mgp = c(2.3, 0.7, 0), las = 1)
    par(cex = 0.85)
    for (S in names(RES)) {
      r <- RES[[S]]
      ev <- r$ds$death == 1
      a1 <- survfit(Surv(r$cs_sb, ev) ~ 1)
      a2 <- survfit(Surv(r$cs_ph, ev) ~ 1)
      lim <- c(0, quantile(c(r$cs_sb, r$cs_ph), 0.995))
      plot(a1$time, a1$cumhaz, type = "s", col = "darkgreen", lwd = 1.4, xlim = lim, ylim = lim,
           xlab = "Cox-Snell residual", ylab = if (S == "AA_R") "Nelson-Aalen cumulative hazard" else "",
           main = c(AA_R = "AA, Regional", WA_R = "WA, Regional")[S])
      lines(a2$time, a2$cumhaz, type = "s", col = "red", lwd = 1.4)
      abline(0, 1, lty = 2)
      if (S == "AA_R")
        legend("topleft", c("Spatial SBART", "PH-frailty"), col = c("darkgreen", "red"), lwd = 1.4, bty = "n")
    }
    invisible(dev.off())
  }
  lab <- c(Age = "Age (years)", HR = "HR status", Grade2 = "Grade 2", Grade3 = "Grade 3",
           TD = "Short TD", BD = "Short BD")
  for (S in names(RES)) {
    r <- RES[[S]]
    yl <- range(-3, 3, r$dv_sb, r$dv_ph)
    fig(sprintf("FigureS%d_deviance_%s.pdf", if (S == "AA_R") 7 else 8, S), 9, 6.6)
    par(mfrow = c(3, 4), mar = c(3.4, 3.6, 2, 0.5), mgp = c(2.1, 0.6, 0), las = 1)
    par(cex = 0.7)
    xv <- with(r$ds, list(Age = Age, HR = as.integer(HR == 1), Grade2 = as.integer(Tg == 2),
                          Grade3 = as.integer(Tg == 3), TD = as.integer(TX == 1), BD = as.integer(BX == 1)))
    for (v in names(lab)) for (m in c("sb", "ph")) {
      dv <- r[[paste0("dv_", m)]]
      ttl <- c(sb = "Spatial SBART", ph = "PH-frailty")[m]
      if (v == "Age") {
        plot(xv$Age, dv, pch = 16, cex = 0.35, ylim = yl, xlab = lab[v], ylab = "Deviance residual", main = ttl)
      } else {
        boxplot(dv ~ factor(xv[[v]], levels = 0:1), ylim = yl, xlab = lab[v], ylab = "Deviance residual",
                main = ttl, pch = 16, cex = 0.35)
      }
      abline(h = c(-2, 2), col = "red", lty = 2)
    }
    invisible(dev.off())
  }

  if ("AA_R" %in% have) {
    P <- R$AA_R
    ct <- CT$AA_R
    ## ---- Figures S9-S10: maps of the county DMS(10), AA Regional
    m <- maps::map("county", "florida", fill = TRUE, plot = FALSE)
    key <- function(s) gsub("[^a-z]", "", tolower(s))
    mk <- key(sub(":.*", "", sub("^florida,", "", m$names)))
    pc <- match(mk, key(COUNTY))
    pc[is.na(pc)] <- match(mk[is.na(pc)], key(sub("Miami-Dade", "Dade", COUNTY)))
    stopifnot(!anyNA(pc))
    draw_map <- function(file, v, brk, labs) {
      val <- rep(NA_real_, N)
      val[ct$county] <- v
      pal <- grey.colors(length(brk) - 1, start = 0.82, end = 0.25)
      col <- pal[cut(val[pc], brk, include.lowest = TRUE)]
      col[is.na(col)] <- "white"
      fig(file, 6.2, 5.2)
      par(mar = c(0.5, 0.5, 0.5, 0.5))
      maps::map("county", "florida", fill = TRUE, col = col, lwd = 0.4)
      legend("bottomleft", c(labs, "no AA Regional patients"), fill = c(pal, "white"), title = "DMS(10), years",
             bty = "n", cex = 0.8)
      invisible(dev.off()) }
    ## SMU: the contrast is exactly 0 in the county with the largest M_hat, which gets a class of its own
    z <- abs(ct$SMU10) <= 1e-9
    brk <- pretty(range(ct$SMU10), 5)
    if (any(z)) brk <- sort(c(brk[abs(brk) > 1e-9], -1e-9, 1e-9))
    num <- function(x) ifelse(abs(x) <= 1e-9, "0", as.character(signif(x, 6)))
    zero_class <- abs(head(brk, -1) + 1e-9) < 1e-12 & abs(brk[-1] - 1e-9) < 1e-12
    labs <- ifelse(zero_class, sprintf("0 (%s)", paste(ct$name[z], collapse = ", ")),
                   paste(num(head(brk, -1)), "to", num(brk[-1])))
    draw_map("FigureS9_map_dms10_smu.pdf", ct$SMU10, brk, labs)
    bt <- pretty(range(ct$TD10), 5)
    draw_map("FigureS10_map_dms10_td.pdf", ct$TD10, bt, sprintf("%.3f to %.3f", head(bt, -1), bt[-1]))

    ## ---- Figures S11-S12: county DMS(10) for TD, AA Regional
    plot_td_by_n("FigureS11_dms10_td.pdf", ct, 0.8, range(ct$TD10_lo, ct$TD10_hi), "DMS(10)")
    sel <- ct[ct$n_i < 5 | (ct$n_i > 50 & ct$n_i < 100), ]
    sel <- sel[order(sel$n_i, sel$county), ]
    fig("FigureS12_dms10_td_ci.pdf", 6.5, 4.1)
    par(mar = c(6.5, 4.4, 0.8, 0.8), las = 1)
    xs <- seq_len(nrow(sel))
    plot(xs, sel$TD10, ylim = range(0, sel$TD10_lo, sel$TD10_hi), pch = 16, xaxt = "n", xlab = "", ylab = "DMS(10)")
    arrows(xs, sel$TD10_lo, xs, sel$TD10_hi, angle = 90, code = 3, length = 0.03)
    abline(h = 0, lty = 3)
    axis(1, at = xs, labels = sprintf("%s (%d)", sel$name, sel$n_i), las = 2, cex.axis = 0.7)
    invisible(dev.off())

    ## ---- Figures S13-S15: Step 1 for M (AA BRFSS), CAR (runs) vs an independent model
    ## logit(M_i) = mu + e_i, e_i iid N(0, s0), fitted here like Step 1
    X1 <- P$X[[1]]
    bf <- read_brfss("mammo_county_16A.xlsx", N)
    m0 <- bf$m_0
    n0 <- bf$n_0i
    stopifnot(isTRUE(all.equal(n0, X1$n_0i)))
    n1 <- X1$settings[["n1"]]
    b1 <- X1$settings[["b1"]]
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260816)
    u <- qlogis((m0 + 0.5)/(n0 + 1))
    s0 <- 1
    mu <- mean(u)
    keep <- matrix(NA_real_, n1 - b1, N)
    for (it in 1:n1) {
      for (i in 1:N) {
        ui <- rnorm(1, u[i], 0.3)
        la <- dbinom(m0[i], n0[i], plogis(ui), log = TRUE) - dbinom(m0[i], n0[i], plogis(u[i]), log = TRUE) -
              0.5*((ui - mu)^2 - (u[i] - mu)^2)/s0
        if (log(runif(1)) < la) u[i] <- ui
      }
      mu <- rnorm(1, mean(u), sqrt(s0/N))
      s0 <- 1/rgamma(1, 1 + N/2, 1 + sum((u - mu)^2)/2)
      if (it > b1) keep[it - b1, ] <- u
    }
    iid <- list(M_hat = colMeans(plogis(keep)), var = apply(keep, 2, var))
    car <- list(M_hat = X1$M_hat, var = X1$var_logitM, q = X1$M_post_q)
    fig("FigureS13_var_logitM.pdf", 5.6, 4.3)
    par(mar = c(4.2, 4.6, 0.8, 0.8), las = 1)
    plot(log(n0), car$var, pch = 1, ylim = range(0, car$var, iid$var), xlab = expression(log(n[0 * i])),
         ylab = expression(widehat(Var) * "{logit(" * M[i] * ") | " * D[0] * "}"))
    points(log(n0), iid$var, pch = 2, col = "red")
    legend("bottomleft", c("Spatial CAR", "Independent"), pch = 1:2, col = c("black", "red"), bty = "n")
    invisible(dev.off())
    fig("FigureS14_Mhat_comparison.pdf", 4.8, 4.6)
    par(mar = c(4.2, 4.4, 0.8, 0.8), las = 1)
    lim <- range(car$M_hat, iid$M_hat)
    plot(car$M_hat, iid$M_hat, xlim = lim, ylim = lim, xlab = expression(hat(M)[i] * ", spatial CAR model"),
         ylab = expression(hat(M)[i] * ", independent model"))
    abline(0, 1, lty = 2)
    invisible(dev.off())
    fig("FigureS15_Mhat_ci.pdf", 5.6, 4.3)
    par(mar = c(4.2, 4.4, 0.8, 0.8), las = 1)
    plot(log(n0), car$M_hat, pch = 16, cex = 0.7, ylim = range(car$q), xlab = expression(log(n[0 * i])),
         ylab = expression(hat(M)[i] * " and 95% interval"))
    arrows(log(n0), car$q[1, ], log(n0), car$q[3, ], angle = 90, code = 3, length = 0.02, col = "grey40")
    invisible(dev.off())

    ## ---- Figure S16: frailty W_i of 3 counties, AA Regional
    n_i <- tabulate(P$county, nbins = N)
    Wm <- colSums(P$w * P$W_draws)/sum(P$w)
    Wq <- sapply(1:N, function(i) wq(P$W_draws[, i], P$w, c(0.025, 0.975)))
    pick <- c(9, 17, 6)                              # Citrus, Escambia, Broward: small, moderate and large n_i
    pick <- pick[n_i[pick] > 0]
    fig("FigureS16_W_ci.pdf", 5.6, 4.3)
    par(mar = c(4.2, 4.4, 0.8, 0.8), las = 1)
    plot(n_i[pick], Wm[pick], pch = 16, log = "x", xlim = range(n_i[pick]) * c(0.7, 1.4),
         ylim = range(Wq[, pick]) * c(1, 1.08), xlab = expression(n[i] * " (log scale)"),
         ylab = expression(hat(W)[i] * " and 95% interval"))
    arrows(n_i[pick], Wq[1, pick], n_i[pick], Wq[2, pick], angle = 90, code = 3, length = 0.05)
    text(n_i[pick], Wq[2, pick], COUNTY[pick], pos = 3, cex = 0.8)
    invisible(dev.off())

    ## ---- Figure S17: Step-2 traces of the AA Regional chains
    XD <- P$X
    bi <- XD[[1]]$settings[["burn_in"]]
    ni <- XD[[1]]$settings[["num_iter"]]
    THIN <- XD[[1]]$settings[["THIN"]]
    same_length <- function(x) x$settings[["burn_in"]] == bi && x$settings[["num_iter"]] == ni &&
                                 x$settings[["THIN"]] == THIN
    stopifnot("Figure S17 needs AA Regional chains of the same length" = all(sapply(XD, same_length)))
    post <- (bi + 1):ni
    th <- seq(1, length(post), by = THIN)
    cols <- rep(c("blue", "red", "darkgreen"), length.out = length(XD))
    panel <- function(x, L, ylim, ylab, main, tag) {
      plot(x, L[[1]], type = "l", col = cols[1], lwd = 0.5, ylim = ylim, xlab = "Iteration (post burn-in)",
           ylab = ylab, main = main)
      for (k in seq_along(L)[-1]) lines(x, L[[k]], col = cols[k], lwd = 0.5)
      mtext(sprintf("(%s)", tag), side = 3, adj = 0, line = 0.3, font = 2) }
    trace_panel <- function(f, ylab, main, tag) {
      L <- lapply(XD, function(x) f(x)[post])
      panel(post[th] - bi, lapply(L, `[`, th), range(unlist(L)), ylab, main, tag) }
    TRC <- c(9, 35, 37)                                  # Citrus, Lake, Leon
    fig("FigureS17_traces.pdf", 9, 5.2)
    par(mfrow = c(2, 3), mar = c(3.6, 4.2, 2, 0.6), mgp = c(2.4, 0.7, 0), las = 1)
    par(cex = 0.75)
    trace_panel(function(x) log(x$s1_trace), expression(log(sigma[1]^2)), "", "a")
    trace_panel(function(x) x$rho1_trace, expression(rho[1]), "", "b")
    t3 <- which.min(abs(XD[[1]]$times - 3))
    S3 <- lapply(XD, function(x) x$S_rep_draws[t3, ])                   # one draw per THIN iterations
    panel(1 + (seq_along(S3[[1]]) - 1) * THIN, S3, range(unlist(S3)), "S(3)", "Broward: 3-year survival", "c")
    for (j in 1:3)
      trace_panel(function(x) x$R_trace[, TRC[j]], bquote(log(W[.(TRC[j])])), COUNTY[TRC[j]], letters[3 + j])
    invisible(dev.off())
  }

  ## ==== paper_numbers.txt: the numbers quoted in the text
  txt <- character(0)
  say <- function(...) txt <<- c(txt, sprintf(...))
  lst <- function(x, sep = "; ") paste(names(x), x, collapse = sep)          # "name value; name value; ..."
  per <- function(f) setNames(sapply(have, f), have)                         # one value per stratum
  chains <- function(s) {
    X <- R[[s]]$X
    setting <- function(k) sapply(X, function(x) x$settings[[k]])
    paste(sprintf("%s (%d iterations, burn-in %d, %d draws)", names(X), setting("num_iter"), setting("burn_in"),
                  sapply(X, `[[`, "n_keep")), collapse = ", ") }
  pval <- function(chisq) pchisq(chisq, 1, lower.tail = FALSE)
  say("Numbers of the data analysis: Step-3 weighted posterior summaries, chains pooled; DMS in years;")
  say("TD and BD contrasts are short minus long delay; intervals are 95%% HPD unless stated otherwise.")
  say("Chains pooled:")
  for (s in have) say("  %s: %s", s, chains(s))

  say("\n== Main text, Section 5 ==")
  say("Counties: %d", N)
  say("Patients (events) per stratum: %s", lst(per(function(s) sprintf("%d (%d)", R[[s]]$n, sum(R[[s]]$delta)))))
  say("Stratified log-rank test AA vs WA (strata = stage): chi-square %.1f on 1 df, p = %.2g",
      LR_strat, pval(LR_strat))
  say("Table 1, VIP (also table1_vip.csv):")
  for (k in names(t1))
    say("  %-16s %s", sub(" \\(\\$t\\$\\)", " t", t1[k]), lst(setNames(cell(VIP[k, ]), STRATA), "  "))
  say("Covariates ranked by VIP (the time t has VIP %s):", lst(setNames(f3(VIP["time", have]), have), ", "))
  for (s in have) {
    v <- sort(VIP[1:6, s], decreasing = TRUE)
    say("  %s: %s", s, paste(names(v), f3(v), collapse = " > "))
  }
  F3m <- sapply(vals, function(x) sapply(x, function(v) if (length(v)) median(v) else NA))
  say(paste("Figure 3, median over the counties, AA / WA: long BD %.3f / %.3f; No SMU (1 - M_hat) %.3f / %.3f;",
            "long TD %.3f / %.3f"), F3m[1, 1], F3m[1, 2], F3m[2, 1], F3m[2, 2], F3m[3, 1], F3m[3, 2])
  for (r in c("AA", "WA")) {
    M <- Mhat_race(r)
    if (length(M))
      say("M_hat, %s: Broward %.3f; median of the %d counties %.3f; range %.3f to %.3f (largest: %s)",
          r, M[BROWARD], N, median(M), min(M), max(M), COUNTY[which.max(M)])
  }
  say("Figure 4, S(t) of the representative patient in Broward, posterior mean (95%% HPD):")
  for (s in have) {
    P <- R[[s]]
    St <- sapply(c(1, 2, 3, 5, 10), function(t)
      sprintf("S(%d) %s", t, est(P$S_rep_draws[which.min(abs(P$times - t)), ], P$w)))
    say("  %s: %s", s, paste(St, collapse = "; "))
  }
  if (all(c("AA_L", "WA_L") %in% have))
    say("Figure 4, Local stage: largest |S_AA(t) - S_WA(t)| over the grid %.3f", max(abs(S_mean$AA_L - S_mean$WA_L)))
  say("Table 2, DMS(10) (95%% HPD), representative AA patient in Broward (also table2_dms10_broward.csv):")
  for (cn in c("SMU", "BD", "TD")) for (g in 1:3) {
    x <- T2[T2$contrast == cn & T2$grade == g, ]
    cells <- ifelse(is.na(x$DMS10), "--", sprintf("%s (%s, %s)", f3(x$DMS10), f3(x$HPD_lo), f3(x$HPD_hi)))
    say("  %-3s grade %d: %s", cn, g, lst(setNames(cells, x$stage)))
  }
  if ("AA_L" %in% have) {
    P <- R$AA_L
    ct <- CT$AA_L
    say("TD, AA Local, Broward: largest DMS(10) over the 3 grades %s (the representative patient has grade %d)",
        f3(max(T2$DMS10[T2$stratum == "AA_L" & T2$contrast == "TD"])), rep_grade(P))
    for (g in 1:3) {
      v <- P$DMS10[, sprintf("TD_g%d", g)]
      say("  grade %d: posterior probability that DMS(10) is > 0: %.3f, = 0: %.3f, < 0: %.3f", g,
          wm(v > 1e-12, P$w), wm(abs(v) <= 1e-12, P$w), wm(v < -1e-12, P$w))
    }
    say(paste("Figure 5 (AA Local, TD): %d counties with patients, n_i from %d to %d",
              "(%d counties with n_i <= 2)"),
        nrow(ct), min(ct$n_i), max(ct$n_i), sum(ct$n_i <= 2))
    say("  DMS(10) estimates %s to %s; lower HPD limits %s to %s; upper HPD limits %s to %s",
        f3(min(ct$TD10)), f3(max(ct$TD10)), f3(min(ct$TD10_lo)), f3(max(ct$TD10_lo)),
        f3(min(ct$TD10_hi)), f3(max(ct$TD10_hi)))
    p_neg <- sapply(seq_len(nrow(ct)), function(i) wm(P$cty_dms[i, "TD10", ] < -1e-12, P$w))
    say("  largest posterior probability of a negative DMS(10) over the counties %.3f", max(p_neg))
  }
  say("PH-frailty hazard ratio (95%% CI, p), short vs long TD:")
  for (s in have) say("  %s: %s", s, ph_hr(s, "I(TX == 1)TRUE"))
  zp <- per(function(s) PH[[s]]$zph["I(TX == 1)", "p"])
  say(paste("PH-frailty, Grambsch-Therneau test of proportional hazards for TD (cox.zph, scaled Schoenfeld residuals,",
            "km transform), p: %s;",
            "largest p = %.3g"), lst(setNames(sprintf("%.3g", zp), names(zp))), max(zp))
  n_excl <- function(s) sum(CT[[s]]$SMU10_lo > 0 | CT[[s]]$SMU10_hi < 0)
  say("Counties whose 95%% HPD interval of the SMU DMS(10) excludes 0: %s",
      lst(per(function(s) sprintf("%d of %d", n_excl(s), nrow(CT[[s]])))))

  say("\n== Supplement, Section 2 ==")
  say("Log-rank test AA vs WA within each stage: %s",
      lst(setNames(sprintf("chi-square %.1f, p = %.2g", unlist(LR), pval(unlist(LR))), names(LR))))
  out2 <- function(s)
    sprintf("%.1f%% vs %.1f%%", 100*mean(abs(R[[s]]$deviance_res) > 2), 100*mean(abs(PH[[s]]$dev) > 2))
  say("Deviance residuals outside (-2, 2), spatial SBART vs PH-frailty: %s", lst(per(out2)))
  for (S in names(RES)) {
    r <- RES[[S]]
    ev <- r$ds$death == 1
    nad <- function(cs) {
      f <- survfit(Surv(cs, ev) ~ 1)
      g <- f$time <= quantile(cs, 0.99)
      max(abs(f$cumhaz[g] - f$time[g])) }
    say(paste("Figure S%d (%s): outside (-2, 2) %.1f%% (SBART) vs %.1f%% (PH-frailty);",
              "above +2 / below -2: SBART %d / %d, PH-frailty %d / %d"),
        if (S == "AA_R") 7 else 8, S, 100*mean(abs(r$dv_sb) > 2), 100*mean(abs(r$dv_ph) > 2),
        sum(r$dv_sb > 2), sum(r$dv_sb < -2), sum(r$dv_ph > 2), sum(r$dv_ph < -2))
    say(paste("Figure S6 (%s): largest |Nelson-Aalen - 45-degree line| up to the 99th percentile of the",
              "residuals: SBART %.3f, PH-frailty %.3f"), S, nad(r$cs_sb), nad(r$cs_ph))
  }
  if ("AA_R" %in% have) {
    ct <- CT$AA_R
    wd <- ct$TD10_hi - ct$TD10_lo
    a <- pick[1]
    b <- pick[length(pick)]
    say("Figure S9 (AA Regional, SMU): county DMS(10) from %s to %s; exactly 0 in %s",
        f3(min(ct$SMU10)), f3(max(ct$SMU10)), paste(ct$name[abs(ct$SMU10) <= 1e-9], collapse = ", "))
    say(paste("Figures S10-S11 (AA Regional, TD): %d counties; DMS(10) estimates from %s to %s;",
              "HPD interval excludes 0 in %d counties"),
        nrow(ct), f3(min(ct$TD10)), f3(max(ct$TD10)), sum(ct$TD10_lo > 0 | ct$TD10_hi < 0))
    grp <- list(`n_i < 5` = ct$n_i < 5, `50 < n_i < 100` = ct$n_i > 50 & ct$n_i < 100,
                `n_i >= 100` = ct$n_i >= 100)
    say("Figure S12: median width of the HPD interval of the TD DMS(10): %s",
        lst(sapply(grp, function(j) sprintf("(%d counties) %.3f", sum(j), median(wd[j])))))
    W3 <- sprintf("n_i %s: %.2f (%.2f, %.2f; %.2f)", ifelse(n_i[pick] < 10, "< 10", paste("=", n_i[pick])), Wm[pick], Wq[1, pick], Wq[2, pick],
                  Wq[2, pick] - Wq[1, pick])
    say("Figure S16 (AA Regional), W_i posterior mean (95%% interval; width): %s", lst(setNames(W3, COUNTY[pick])))
    if (length(pick) >= 2)
      say("Figure S16: width ratio %s / %s %.2f, versus %.2f = sqrt of the ratio of their n_i",
          COUNTY[a], COUNTY[b], diff(Wq[, a])/diff(Wq[, b]), sqrt(n_i[b]/n_i[a]))
    sm <- n0 <= quantile(n0, 0.25)
    lq <- n0 >= quantile(n0, 0.75)
    ciw <- car$q[3, ] - car$q[1, ]
    say("Figure S13 (AA BRFSS): n_0i from %s to %d; Var{logit(M_i) | D_0} under the CAR model from %.3f to %.3f",
        if (min(n0) < 10) "< 10" else min(n0), max(n0), min(car$var), max(car$var))
    say(paste("Figure S13: median Var{logit(M_i) | D_0}, lowest / highest quartile of n_0i:",
              "CAR %.3f / %.3f, independent model %.3f / %.3f"),
        median(car$var[sm]), median(car$var[lq]), median(iid$var[sm]), median(iid$var[lq]))
    say(paste("Figures S13-S14: posterior variance larger under the independent model in %d of %d counties;",
              "largest |M_hat(independent) - M_hat(CAR)| %.3f"),
        sum(iid$var > car$var), N, max(abs(iid$M_hat - car$M_hat)))
    say(paste("Figure S15: median width of the 95%% interval of M_i, lowest / highest quartile of n_0i:",
              "%.3f / %.3f; Spearman correlation of M_hat_i and log(n_0i) %.2f"),
        median(ciw[sm]), median(ciw[lq]), cor(log(n0), car$M_hat, method = "spearman"))
    say("Figure S17 (AA Regional): %d chains (seeds %s), %d iterations, burn-in %d",
        length(XD), paste(sapply(XD, function(x) x$settings[["seed"]]), collapse = ", "), ni, bi)
  }
  writeLines(txt, file.path(figdir, "paper_numbers.txt"))
  cat("wrote", length(list.files(figdir)), "files to", figdir, "\n")
}

################################################################################
# Command line
################################################################################
args <- commandArgs(trailingOnly = TRUE)
cmd <- if (length(args)) args[1] else ""
if (cmd == "fit" && length(args) >= 4) {
  fit_stratum(args[-1])
} else if (cmd == "report" && length(args) == 3) {
  make_report(args[2], args[3])
} else {
  cat("Usage:\n",
      "  Rscript fcr_sbart_data_analysis.R fit <RACE> <STAGE> <outdir> [num_iter burn_in num_tree use_step3\n",
      "          THIN NQ n1 b1 spatial county t_horizon seed]\n",
      "  Rscript fcr_sbart_data_analysis.R report <runs_dir> <figdir>\n", sep = "")
  if (!interactive()) quit(status = 1)
}
