
################################################################################
# FCR data analysis -- ALL THREE STEPS of the algorithm, one stratum per run.
#
#   Step 1  CAR model for the county-level unobservable M from the BRFSS data
#           D_0 alone; MH for logit(M_i), Gibbs for sigma_0^2, MH for rho_0.
#           Returns the posterior SAMPLE M_post and the Bayes estimate
#           M_hat = E[M | D_0].
#   Step 2  SBART survival model on D_1 with Poisson-thinning data augmentation,
#           M_hat plugged in.  Gibbs for lambda_0 and sigma_1^2, MH for rho_1,
#           component-wise MH for R = log W.  This targets the APPROXIMATE
#           posterior with M fixed at M_hat.
#   Step 3  Importance reweighting.  Each retained Step-2 draw is paired with a
#           fresh draw M^(b) from p_0(M | D_0) and given the weight
#             log omega* = sum_ij [ delta_ij log{Phi(b(y_ij,M_i))/Phi(b(y_ij,Mhat_i))}
#                          - lambda_0 W_i int_0^{y_ij}{Phi(b(s,M_i))-Phi(b(s,Mhat_i))} ds ]
#           which converts the Step-2 sample into a sample from the exact joint
#           posterior.  EVERY posterior summary below is a self-normalised
#           weighted average over these draws, and the importance-sampling ESS
#           is reported so the reader can judge them.
#
# Outputs, all computed under the 3 step posterior 
#   * Variable Inclusion Proportions (Table 1)
#   * survival curves S(t) for a representative patient (Figure 4)
#   * DMS(10) and DMS(5) with 95% HPD by tumour grade (Table 2, Table S2)
#   * DMS(10) and DMS(5) for the representative patient of EVERY county
#     (Figure 5, Figures S9-S12)
#   * Cox-Snell and deviance residuals (Figures S6-S8)
#   * M_hat, Var{logit(M_i) | D_0} and 95% intervals of M_i (Figures S13, S15)
#   * draws and 95% intervals of the frailties W_i (Figure S16)
#   * full MCMC traces of lambda_0, sigma_1^2, rho_1 and log W (Figure S17)
# Two files are written to outdir: fcr_<race>_<stage>.rds (summaries) and
# fcr_<race>_<stage>_details.rds (everything make_figures_tables.R needs).
#
# Data files read from the working directory (not distributed):
#   Surv_data2.csv, mammo_county_16A.xlsx, mammo_county_16W.xlsx, W.mat.csv
#

#

################################################################################
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)
race   <- if (length(args) >= 1) toupper(args[1]) else "AA"
stage  <- if (length(args) >= 2) toupper(args[2]) else "D"
outdir <- if (length(args) >= 3) args[3] else sprintf("fcr_%s_%s", race, stage)
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

suppressMessages({library(SoftBart); library(MASS); library(truncnorm); library(readxl)})
logf <- file.path(outdir, "analysis.log")
lg <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), " [", race, "-", stage, "] ", ...)
                      cat(m, "\n", file = logf, append = TRUE); cat(m, "\n"); flush.console() }
t_start <- Sys.time()
RNGkind("L'Ecuyer-CMRG")
set.seed(20260816)   # Step 1 always uses this stream, so every chain shares the same M_hat

################################################################################
# 0. Data
################################################################################
A <- as.matrix(read.csv("W.mat.csv", row.names = 1))
N <- nrow(A); D <- diag(rowSums(A))
stopifnot(county_oi %in% seq_len(N))
## the adjacency must be in the county order of the data (FIPS order): check some
## known Florida neighbours (Broward 6 - Palm Beach 50, Miami-Dade 13 - Monroe 44,
## Escambia 17 - Santa Rosa 57, Hillsborough 29 - Pinellas 52, Leon 37 - Gadsden 20)
adj_ok <- N == 67 && all(A[cbind(c(6, 13, 17, 29, 37), c(50, 44, 57, 52, 20))] == 1)

surv_raw <- read.csv("Surv_data2.csv")
names(surv_raw) <- c("id","time_days","county","TX_Delay","death",
                     "BX_Delay","Age","HR_p","Tgrade","Race","SMUprob","Stage")
race_code  <- if (race == "AA") 2 else 1          # Race: 2 = African American
stage_code <- switch(stage, L = 1, R = 2, D = 3)
st <- surv_raw[surv_raw$Race == race_code & surv_raw$Stage == stage_code, ]

## Follow-up in years.  With an administrative horizon T_ADMIN, a subject still
## at risk at the horizon is censored there (delta set to 0).
T_ADMIN <- Inf                                     # no administrative cap by default
yr <- st$time_days / 365.25
delta <- as.integer(st$death == 1)
if (is.finite(T_ADMIN)) { delta[yr > T_ADMIN] <- 0L; yr <- pmin(yr, T_ADMIN) }
yr <- pmax(yr, 1/365.25)                           # a same-day event is 1 day, not 0

sd_df <- data.frame(
  county = st$county, time = yr, delta = delta,
  Age    = st$Age,
  HR     = as.integer(st$HR_p == 1),
  Grade2 = as.integer(st$Tgrade == 2),
  Grade3 = as.integer(st$Tgrade == 3),
  TD     = as.integer(st$TX_Delay == 3),
  BD     = as.integer(st$BX_Delay == 3))
sd_df <- sd_df[!is.na(sd_df$county) & sd_df$county >= 1 & sd_df$county <= N &
               !is.na(sd_df$time) & sd_df$time > 0 & !is.na(sd_df$delta), ]
nt <- nrow(sd_df); cty <- sd_df$county

cov_names <- c("Age","HR","Grade2","Grade3","TD","BD")
Xraw <- as.matrix(sd_df[, cov_names])
## scale continuous covariates to (0,1) as Section 2.2 requires; binary stay 0/1
X <- apply(Xraw, 2, function(x) { r <- range(x, na.rm = TRUE)
  if (diff(r) < 1e-10) rep(0.5, length(x)) else (x - r[1])/diff(r) })
p <- ncol(X)

## BRFSS D_0 for this race.  m_0i = number reporting mammography use, n_0i the
## county total.  A county absent from the file gets n_0i = 0 and contributes NO
## binomial term -- the CAR prior interpolates it.  (Coding it as 0 successes
## out of 1 would inject a false observation.)
bf_file <- if (race == "AA") "mammo_county_16A.xlsx" else "mammo_county_16W.xlsx"
bf <- suppressMessages(read_excel(bf_file, skip = 4, col_names = FALSE))
bcode <- suppressWarnings(as.numeric(bf[[1]]))
## county code = (FIPS + 1)/2; Miami-Dade is FIPS 086 (formerly 025, i.e. county 13)
bcty <- ifelse(!is.na(bcode) & bcode == 86, 13, (bcode + 1)/2)
b_no  <- suppressWarnings(as.numeric(bf[[2]]))     # column "0" : no mammography
b_yes <- suppressWarnings(as.numeric(bf[[3]]))     # column "1" : mammography
ok <- !is.na(bcty) & bcty >= 1 & bcty <= N
m_0 <- rep(0, N); n_0i <- rep(0, N)
m_0[round(bcty[ok])]  <- b_yes[ok]
n_0i[round(bcty[ok])] <- b_no[ok] + b_yes[ok]
has0 <- n_0i > 0
lg(sprintf("stratum n=%d events=%d (%.1f%% censored) | counties with BRFSS data: %d/%d",
           nt, sum(sd_df$delta), 100*mean(sd_df$delta == 0), sum(has0), N))
if (!adj_ok) lg("WARNING: W.mat.csv does not match the county order of the data (known Florida neighbours are not adjacent); the CAR models will use wrong neighbours")
lg(sprintf("settings: K=%d iter=%d burn=%d THIN=%d NQ=%d step3=%d spatial=%d",
           num_tree, num_iter, burn_in, THIN, NQ, use_step3, spatial))

## rho support from the eigenvalues of the normalised adjacency D^{-1/2} A D^{-1/2},
## capped at (-1,1).
dv <- rowSums(A); An <- diag(1/sqrt(dv)) %*% A %*% diag(1/sqrt(dv))
ev_n <- eigen((An + t(An))/2, only.values = TRUE)$values
rho_bounds <- c(max(-1, 1/min(ev_n)), min(1, 1/max(ev_n)))

## log|D - rho A| = sum(log a_i+) + sum_k log(1 - rho e_k), e_k the eigenvalues of
## D^{-1/2} A D^{-1/2}
car_ld <- function(x, rho, s2)
  0.5*(sum(log(dv)) + sum(log(pmax(1 - rho*ev_n, 1e-300))) - N*log(s2)) -
  0.5*as.numeric(t(x) %*% (D - rho*A) %*% x)/s2
## truncated uniform random walk for rho: proposal and window width (Hastings factor)
rho_prop <- function(r, h) runif(1, max(rho_bounds[1], r - h), min(rho_bounds[2], r + h))
rho_win  <- function(r, h) min(rho_bounds[2], r + h) - max(rho_bounds[1], r - h)
iid_ld <- function(x, rho, s2) sum(dnorm(x, 0, sqrt(s2), log = TRUE))
frail_ld <- if (spatial == 1) car_ld else iid_ld     # spatial = 0 is the competitor

################################################################################
# 1. STEP 1 -- CAR model for M from D_0 alone
################################################################################
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

################################################################################
# 2. STEP 2 -- SBART survival with Poisson-thinning data augmentation
################################################################################
lg("=== STEP 2: SBART survival ===")
max_time <- max(sd_df$time)
Xd <- cbind(sd_df$time/max_time, M_hat[cty], X)
colnames(Xd) <- c("time", "M", cov_names)
Yinit <- ifelse(sd_df$delta == 1, 0.5, -0.5) + rnorm(nt, 0, 0.1)
## sigma fixed at 1, as the probit data augmentation requires. k = 1/3 sets the scale
## 3/(2 sqrt(K)) of the leaf standard deviation sigma_mu; the other SBART hyperparameters
## (half-Cauchy update of sigma_mu, splitting proportions, bandwidths) are the SoftBart
## defaults of Linero and Yang (2018), including their warm-up of num_burn = burn_in.
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
cty_idx <- lapply(1:N, function(c) which(cty == c))

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
## covariate values equal to the sample medians: the median for age, the lower
## median (an observed value) for the binary covariates, and the lower median
## tumour grade, coded by the two grade indicators
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

## ===== Additional prediction scenarios for the tables and figures (do not affect the MCMC) =====
## (a) Broward (county_oi) grade-specific contrasts, for Table 2 / Table S2
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
    ## the exponent of W_i is delta_i+ + m_i+ : every rejected point also carries
    ## a factor lambda_0 W_i in the complete-data likelihood
    llp <- (d_cnt[c2] + m_cnt[c2])*Rp[c2] - lambda0*exp(Rp[c2])*sy[c2]
    llc <- (d_cnt[c2] + m_cnt[c2])*Rv[c2] - lambda0*exp(Rv[c2])*sy[c2]
    la <- (llp + frail_ld(Rp, rho1, s1)) - (llc + frail_ld(Rv, rho1, s1))
    if (is.finite(la) && log(runif(1)) < la) { Rv <- Rp; acc_R[c2] <- acc_R[c2] + 1 }
  }
  Wv <- exp(Rv)          # no deterministic recentring: the CAR here is proper
  keep_lam[iter] <- lambda0; keep_s1[iter] <- s1; keep_rho1[iter] <- rho1
  R_trace[iter, ] <- Rv

  ############################################################################
  # 3. STEP 3 -- importance weight for this draw, and the weighted summaries
  ############################################################################
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
    ## variable inclusion proportions, weighted and post-burn-in only
    ## get_counts() returns a one-column matrix; as.numeric keeps vip_sum a
    ## NAMED vector so the covariates can be addressed by name below
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
    W_draws[n_keep, ] <- Wv; M_draws[n_keep, ] <- Mdraw
    vip_draws[n_keep, ] <- cnt
    if (use_step3 == 1) {
      ## the five county-specific scenarios take the drawn M_i; the "best county"
      ## scenario is a stipulated value and stays as it is
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

################################################################################
# 4. Posterior summaries -- every one of them Step-3 weighted
################################################################################
lwk <- lw_keep[1:n_keep]; Wt <- exp(lwk - max(lwk))
ess_w <- sum(Wt)^2/sum(Wt^2)
w_top <- max(Wt)/sum(Wt); p_zero <- mean(lwk == 0)
if (use_step3 == 1) {
  lg(sprintf("Step 3: %d retained draws, ESS = %.1f (%.0f%%), largest normalised weight %.3f, share of log-weights = 0: %.2f",
             n_keep, ess_w, 100*ess_w/n_keep, w_top, p_zero))
} else lg(sprintf("Step 3 OFF (2-step method, unit weights): %d retained draws", n_keep))
if (ess_w < 0.1*n_keep)
  lg("WARNING: importance weights are close to degenerate; treat interval estimates with caution")

S_rep   <- S_scen_sum/sumw                     # posterior mean survival, per scenario
S_repM  <- if (use_step3 == 1) S_scenM_sum/sumw else NULL
vip_raw <- vip_sum/sumw
## VIP over the paper's 6 covariates: tumour grade is the two dummies combined
## VIP = posterior (Step-3 weighted) mean of the proportion of splitting rules
## that use each variable; `time` is one of the variables, so the six covariates
## need not sum to 1.
vip6 <- c(Age = unname(vip_raw["Age"]), HR = unname(vip_raw["HR"]),
          Grade = unname(vip_raw["Grade2"] + vip_raw["Grade3"]), SMU = unname(vip_raw["M"]),
          TD = unname(vip_raw["TD"]), BD = unname(vip_raw["BD"]))
csnell <- csnell_sum/sumw
mart   <- sd_df$delta - csnell
devres <- sign(mart) * sqrt(pmax(-2*(mart + sd_df$delta*log(pmax(csnell, 1e-12))), 0))

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


################################################################################
# 5. Detailed outputs for make_figures_tables.R (Step-3 weighted AND unweighted = 2-step)
################################################################################
wq <- function(v, w, p) { o <- order(v); cw <- cumsum(w[o])/sum(w); v[o][which(cw >= p)[1]] }
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
cd <- cty_dms[, , 1:n_keep, drop = FALSE]
cty_tab <- do.call(rbind, lapply(seq_len(nC), function(k) {
  row <- data.frame(county = cty_list[k], n_i = n_i[cty_list[k]], M_hat = M_hat[cty_list[k]], n_0i = n_0i[cty_list[k]])
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
vip_norm <- function(w) { pr <- vd/pmax(rowSums(vd), 1); v <- colSums(w*pr)/sum(w)
  c(Age = unname(v["Age"]), HR = unname(v["HR"]), Grade = unname(v["Grade2"] + v["Grade3"]),
    SMU = unname(v["M"]), TD = unname(v["TD"]), BD = unname(v["BD"]), time = unname(v["time"])) }
saveRDS(list(race = race, stage = stage, n = nt, events = sum(sd_df$delta),
  settings = c(num_iter = num_iter, burn_in = burn_in, num_tree = num_tree, THIN = THIN, NQ = NQ,
               n1 = n1, b1 = b1, spatial = spatial, use_step3 = use_step3, county = county_oi, t_hor = t_hor,
               seed = seed),
  ess_weights = ess_w, w_top = w_top, p_logw_zero = p_zero, n_keep = n_keep, logw = lwk,
  tab2 = tab2, cty_tab = cty_tab, W_tab = W_tab,
  vip6_weighted = vip_norm(Wt), vip6_unweighted = vip_norm(U), vip_draws = vd,
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
