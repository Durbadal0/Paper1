

options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)
rep_id   <- as.integer(args[1])
outdir   <- if (length(args) >= 2) args[2] else "scenA_out"
num_iter <- if (length(args) >= 3) as.integer(args[3]) else 10000
burn_in  <- if (length(args) >= 4) as.integer(args[4]) else 3500
num_tree <- if (length(args) >= 5) as.integer(args[5]) else 200
scenario  <- if (length(args) >= 6) toupper(args[6]) else "A"
use_step3 <- if (length(args) >= 7) as.integer(args[7]) else 0
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({library(SoftBart); library(MASS); library(truncnorm)
  library(survival); library(ranger)})

logf <- file.path(outdir, sprintf("rep%02d.log", rep_id))
lg <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), " [rep ", rep_id, "] ", ...)
cat(m, "\n", file = logf, append = TRUE); cat(m, "\n"); flush.console() }
t_start <- Sys.time()



make_A15 <- function() {
  cn <- c("Escambia","SantaRosa","Okaloosa","Walton","Holmes","Washington","Bay",
          "Jackson","Calhoun","Gulf","Gadsden","Liberty","Franklin","Leon","Wakulla")
  ed <- list(c("Escambia","SantaRosa"), c("SantaRosa","Okaloosa"), c("Okaloosa","Walton"),
             c("Walton","Holmes"), c("Walton","Washington"), c("Walton","Bay"),
             c("Holmes","Washington"), c("Holmes","Jackson"),
             c("Washington","Jackson"), c("Washington","Calhoun"), c("Washington","Bay"),
             c("Bay","Calhoun"), c("Bay","Gulf"),
             c("Jackson","Calhoun"), c("Jackson","Gadsden"), c("Jackson","Liberty"),
             c("Calhoun","Gulf"), c("Calhoun","Liberty"),
             c("Gulf","Liberty"), c("Gulf","Franklin"),
             c("Gadsden","Liberty"), c("Gadsden","Leon"),
             c("Liberty","Franklin"), c("Liberty","Leon"), c("Liberty","Wakulla"),
             c("Franklin","Wakulla"), c("Leon","Wakulla"))
  A <- matrix(0, length(cn), length(cn), dimnames = list(cn, cn))
  for (e in ed) { A[e[1], e[2]] <- 1; A[e[2], e[1]] <- 1 }
  A }

make_A10 <- function() { N <- 10; A <- matrix(0, N, N)
for (e in list(c(1,2),c(2,3),c(3,4),c(4,5),c(4,7),c(5,6),c(5,8),
               c(6,7),c(6,8),c(6,9),c(7,9),c(7,10),c(8,9),c(9,10)))
{ A[e[1],e[2]] <- 1; A[e[2],e[1]] <- 1 }; A }

if (scenario == "D") {                       # Scenario D: N = 10 clusters, 1,000 subjects
  A <- make_A10(); n_i <- c(80,100,120,100,80,120,80,100,120,100)      # sum = 1000
  n_0i <- c(8,12,15,20,10,18,9,25,14,11); t_max <- 8
} else {                                     # Scenarios A, B, C: N = 15, 1,250 subjects
  A <- make_A15(); n_i <- c(60,70,80,100,60,90,70,110,90,60,70,90,100,80,120)
  n_0i <- c(8,10,12,15,9,14,11,20,13,8,10,16,25,12,18); t_max <- 10
}
N <- nrow(A); D <- diag(rowSums(A))
p       <- 15          # observed covariates, only x1,x2,x3 active
sig0    <- 1; sig1 <- 1; rho_true <- 0.5
cens_target <- 0.12    # 10-15% censoring

# rho support: normalized adjacency, capped at (-1,1) [corrected manuscript]
dvec <- rowSums(A); Anorm <- diag(1/sqrt(dvec)) %*% A %*% diag(1/sqrt(dvec))
ev_n <- eigen((Anorm + t(Anorm))/2, only.values = TRUE)$values
rho_bounds <- c(max(-1, 1/min(ev_n)), min(1, 1/max(ev_n)))

# piecewise hazard of Scenario A: I1=(0,2), I2=(2,6), I3=(6,Inf)
g1 <- function(t,M,x1,x2,x3) 0.50*t*x1*x2 + 0.30*t*M*x1 + 0.15*M^2 + 0.40*x2*x3*M
g2 <- function(t,M,x1,x2,x3) 0.35*t*x1*x2*x3*M + 0.25*t*M^2 + 0.20*x1*M + 0.15*M^2
g3 <- function(t,M,x1,x2,x3) 0.30*log(t+1)*x1*x2*x3*M + 0.25*t*x2*x3 + 0.20*M^2
## Scenario D replaces the piecewise exponential hazard by a probit-link hazard
g4 <- function(t,M,x1,x2,x3) 0.5*t*x1*x2 + 0.3*t*M*x1 + 0.15*M^2 +
  0.08*t^2*x1*x2*x3 + 0.3*t*x2*x3*M + 0.4*sqrt(x3)*M*t
haz <- function(t, W, M, x1, x2, x3) {
  if (scenario == "D") return(W * pnorm(g4(t,M,x1,x2,x3)))
  v <- ifelse(t < 2, g1(t,M,x1,x2,x3), ifelse(t < 6, g2(t,M,x1,x2,x3), g3(t,M,x1,x2,x3)))
  if (scenario == "B") { W^(1 + 0.2*t) * exp(v) } else { W * exp(v) } }

# cumulative hazard on a fine grid, vectorised over t (trapezoid)
FINE <- seq(0, 30, length.out = 6001); dFINE <- FINE[2] - FINE[1]
cumH_fine <- function(W, M, x1, x2, x3) {
  h <- haz(FINE, W, M, x1, x2, x3)
  c(0, cumsum((h[-1] + h[-length(h)])/2 * dFINE)) }
S_true_fun <- function(tt, W, M, x1, x2, x3) {
  H <- cumH_fine(W, M, x1, x2, x3); exp(-approx(FINE, H, xout = tt, rule = 2)$y) }


set.seed(20260814 + rep_id)
Qinv <- solve(D - rho_true * A + diag(1e-8, N))
R_true <- mvrnorm(1, rep(0, N), sig1^2 * Qinv); R_true <- R_true - mean(R_true)
W_true <- exp(R_true)
## Scenarios C and D: strong within-cluster association between M_i and W_i
M_true <- if (scenario %in% c("C","D")) {
  plogis(rnorm(N, 0.9 * R_true, sqrt(1 - 0.9^2)))
} else {
  plogis(as.numeric(mvrnorm(1, rep(0, N), sig0^2 * Qinv)))
}
m_0    <- rbinom(N, n_0i, M_true)
Mplug  <- m_0 / n_0i                       # naive plug-in used by RSF / PH

nt <- sum(n_i); cty <- rep(1:N, times = n_i)
X  <- matrix(runif(nt * p), nt, p); colnames(X) <- paste0("x", 1:p)
Tobs <- numeric(nt)
for (k in 1:nt) {
  i <- cty[k]; E <- rexp(1)
  H <- cumH_fine(W_true[i], M_true[i], X[k,1], X[k,2], X[k,3])
  Tobs[k] <- if (max(H) < E) max(FINE) else approx(H, FINE, xout = E, ties = "ordered")$y
}

med   <- median(Tobs)
uexp  <- rexp(nt)                                    # independent of Tobs
gap   <- function(r) mean(Tobs > med + uexp/r) - cens_target
r_cal <- tryCatch(uniroot(gap, c(1e-4, 1e4))$root, error = function(e) 1)
Cens  <- med + uexp/r_cal
delta <- as.integer(Tobs <= Cens)
Tobs  <- pmin(Tobs, Cens)
lg(sprintf("scenario %s | step3=%d | data: n=%d events=%d censored=%.1f%% median T=%.2f", scenario, use_step3,
           nt, sum(delta), 100*mean(delta == 0), med))

sd_df <- data.frame(county = cty, time = Tobs, delta = delta, X)

################################################################################
# 3. Evaluation grids (truth pre-computed)
################################################################################
tg   <- seq(t_max/100, t_max, length.out = 100); nT <- length(tg)  # AMSE times
gv   <- seq(0.05, 0.95, length.out = 8)                            # 8x8x8 grid
xp   <- as.matrix(expand.grid(x1 = gv, x2 = gv, x3 = gv)); nXP <- nrow(xp)
S_true_amse <- array(0, c(N, nXP, nT))
for (i in 1:N) for (r in 1:nXP)
  S_true_amse[i, r, ] <- S_true_fun(tg, W_true[i], M_true[i], xp[r,1], xp[r,2], xp[r,3])

# AES: cluster 1, four (x1, p0) panels, x2 = x3 = 0.5.
# Plotted on (0, 3]: under the Scenario A hazard the true survival is ~0 by t = 3.
t_aes <- 3
pt <- seq(t_aes/100, t_aes, length.out = 100); nPT <- length(pt); dtp <- pt[1]
aes_cases <- list(c(x1=0.3, p0=0.75), c(x1=0.7, p0=0.75),
                  c(x1=0.3, p0=0.50), c(x1=0.7, p0=0.50))
S_true_aes <- t(sapply(aes_cases, function(a)
  S_true_fun(pt, W_true[1], a["p0"], a["x1"], 0.5, 0.5)))

# The AES grid above is also the grid on which every 95% interval is built (see
# section 8).  Version 6's separate 10,000-point (t, x, M_1) grid existed only to
# fill the Table-S1 coverage tally and is gone, together with the truth array it
# fed: version 7 computes no coverage.

################################################################################
# 4. Competing methods: RSF and PH-frailty (plug-in M = m0i/n0i)
################################################################################
rsf_df <- cbind(sd_df[, c("time","delta")], as.data.frame(X), M = Mplug[cty])
rsf <- tryCatch(ranger(Surv(time, delta) ~ ., data = rsf_df, num.trees = 1000,
                       min.node.size = 15, seed = 42 + rep_id), error = function(e) NULL)
rsf_S <- function(xmat, Mv, tt) {              # rows of xmat = covariate vectors
  nd <- as.data.frame(cbind(xmat, M = Mv)); names(nd) <- c(paste0("x",1:p), "M")
  pr <- predict(rsf, data = nd)
  t(apply(pr$survival, 1, function(s) approx(pr$unique.death.times, s, xout = tt, rule = 2)$y)) }

## FIX 3: sparse = FALSE keeps the cluster-level M identified.  With the default
## sparse frailty, coxph aliases M to NA and the old code silently set beta_M = 0.
ph_fml <- as.formula(paste("Surv(time, delta) ~", paste(c(paste0("x",1:p), "M"), collapse = " + "),
                           "+ frailty(county, distribution = 'gaussian', sparse = FALSE)"))
ph_dat <- cbind(sd_df[, c("time","delta","county")], as.data.frame(X), M = Mplug[cty])
ph <- tryCatch(coxph(ph_fml, data = ph_dat), error = function(e) NULL)
if (!is.null(ph)) {
  bph <- coef(ph)
  if (any(is.na(bph))) lg(sprintf("WARNING: aliased PH coefficients: %s",
                                  paste(names(bph)[is.na(bph)], collapse = ",")))
  bph[is.na(bph)] <- 0
  frail  <- ph$frail
  frailv <- if (!is.null(frail) && length(frail) == N) as.numeric(frail) else rep(0, N)
  W_ph   <- exp(frailv)
  ## FIX 4: Breslow baseline WITH the frailties in the risk set:
  ##   H0(t) = sum_{t_k <= t} d_k / sum_{i: y_i >= t_k} exp(x_i'beta + frail_i)
  lp_dat <- as.numeric(as.matrix(ph_dat[, paste0("x", 1:p)]) %*% bph[paste0("x", 1:p)]) +
    bph["M"] * ph_dat$M + frailv[ph_dat$county]
  o   <- order(ph_dat$time)
  tso <- ph_dat$time[o]; dso <- ph_dat$delta[o]; rso <- exp(lp_dat)[o]
  atrisk <- rev(cumsum(rev(rso)))                    # sum of risk scores with y >= t
  et  <- unique(tso[dso == 1])
  dk  <- as.numeric(table(factor(tso[dso == 1], levels = et)))
  den <- atrisk[match(et, tso)]
  H0_step <- cumsum(dk / den)
  H0 <- function(tt) approx(et, H0_step, xout = tt, method = "constant",
                            rule = 2, f = 0, ties = "ordered")$y
  lg(sprintf("PH fitted: beta_M = %.3f | %d event times | H0(max t) = %.3f",
             bph["M"], length(et), H0(max(sd_df$time))))
}

################################################################################
# 5. Step 1: CAR model for M from BRFSS data
################################################################################
car_ld <- function(x, rho, s2) {
  Q <- (D - rho * A); ev <- eigen(Q, only.values = TRUE, symmetric = TRUE)$values
  0.5 * (sum(log(pmax(ev, 1e-12))) - N * log(s2)) - 0.5 * as.numeric(t(x) %*% Q %*% x) / s2 }
n1 <- 12000; b1 <- 6000
logitM <- qlogis((m_0 + 0.5)/(n_0i + 1)); s0 <- 1; rho0 <- 0
Msave <- matrix(NA, n1, N)
for (it in 1:n1) {
  for (i in 1:N) {
    prop <- logitM; prop[i] <- rnorm(1, logitM[i], 0.3)
    la <- (car_ld(prop, rho0, s0) + dbinom(m_0[i], n_0i[i], plogis(prop[i]), log = TRUE)) -
      (car_ld(logitM, rho0, s0) + dbinom(m_0[i], n_0i[i], plogis(logitM[i]), log = TRUE))
    if (is.finite(la) && log(runif(1)) < la) logitM <- prop
  }
  qf <- as.numeric(t(logitM) %*% (D - rho0 * A) %*% logitM)
  s0 <- 1/rgamma(1, 1 + N/2, 1 + qf/2)
  rp <- runif(1, max(rho_bounds[1], rho0 - 0.1), min(rho_bounds[2], rho0 + 0.1))
  if (log(runif(1)) < car_ld(logitM, rp, s0) - car_ld(logitM, rho0, s0)) rho0 <- rp
  Msave[it, ] <- plogis(logitM)
}
M_post <- Msave[(b1+1):n1, ]        # Step-1 posterior sample, used by Step 3
M_hat  <- colMeans(M_post)
lg(sprintf("step1 done: cor(M_true, M_hat) = %.3f  (naive plug-in cor = %.3f)",
           cor(M_true, M_hat), cor(M_true, Mplug)))

################################################################################
# 6. Step 2: SBART survival model
################################################################################
max_time <- max(sd_df$time)
Xd <- cbind(sd_df$time/max_time, M_hat[cty], X)          # (time, M, x1..xp)
colnames(Xd) <- c("time", "M", paste0("x", 1:p))
Yinit <- ifelse(delta == 1, 0.5, -0.5) + rnorm(nt, 0, 0.1)
# sigma = 1 (probit DA) and sigma_mu = 3/(2 sqrt(K))  <=>  k = 1/3
hyp <- Hypers(Xd, Yinit, num_tree = num_tree, sigma_hat = 1, k = 1/3)
forest <- MakeForest(hyp, Opts(num_burn = 0, num_save = 1, update_sigma = FALSE), warn = FALSE)

lambda0 <- 0.1; Rv <- rep(0, N); Wv <- exp(Rv); s1 <- 1; rho1 <- 0
keep_lam <- keep_s1 <- keep_rho1 <- numeric(num_iter); keep_W <- matrix(NA, num_iter, N)
acc_R <- rep(0, N); acc_rho1 <- 0

# fixed prediction blocks (design matrices do not change across iterations)
xfill <- matrix(0.5, 1, p - 3)
tg0 <- c(0, tg); nT0 <- length(tg0)
cumtrap <- function(P, h) {                 # P: rows x (nT+1) values of Phi(b) at tg0
  m <- (P[, -1, drop = FALSE] + P[, -ncol(P), drop = FALSE]) / 2
  t(apply(m, 1, cumsum)) * h }
amse_blocks <- lapply(1:N, function(i) {
  bx <- cbind(xp, xfill[rep(1, nXP), , drop = FALSE])
  cbind(rep(tg0, each = nXP)/max_time, M_hat[i], bx[rep(1:nXP, times = nT0), ]) })
pt0 <- c(0, pt); nPT0 <- length(pt0)
aes_block <- do.call(rbind, lapply(aes_cases, function(a)
  cbind(pt0/max_time, a["p0"], a["x1"], 0.5, 0.5, xfill[rep(1, nPT0), , drop = FALSE])))
# same block with column 2 (M_1) overwritten by the current draw of M_1 from
# p(M | D_0); used for the M-integrated credible band of section 8
aesM_base <- aes_block

THIN <- 3
n_keep_max <- sum(seq_len(num_iter) > burn_in & seq_len(num_iter) %% THIN == 0)
S_amse_sum <- array(0, c(N, nXP, nT)); S_aes_sum <- matrix(0, 4, nPT); n_acc <- 0
## EVERY retained draw of the AES curve is stored with its Step-3 log weight, so
## the credible band and the point estimate are computed from the same sample and
## the same weights.  4 x 100 x ~2200 doubles is a few MB.
S_aes_draws  <- array(NA_real_, c(4, nPT, n_keep_max))
S_aesM_draws <- if (use_step3 == 1) array(NA_real_, c(4, nPT, n_keep_max)) else NULL
lw_keep <- numeric(n_keep_max); n_keep <- 0
## --- Step 3 (version 3 only): importance weights omega* -----------------------
sumw <- 0; maxlw <- -Inf; lw_all <- c(); w_lam <- 0; w_W <- rep(0, N)
NQ <- 50; qw <- (seq_len(NQ) - 0.5)/NQ
q_rows  <- rep(1:nt, times = NQ)
Xq_base <- cbind(as.numeric(outer(sd_df$time, qw))/max_time, 0, X[q_rows, , drop = FALSE])
Xy_base <- cbind(sd_df$time/max_time, 0, X)
rescale_all <- function(f) { S_amse_sum <<- S_amse_sum*f; S_aes_sum <<- S_aes_sum*f
w_lam <<- w_lam*f; w_W <<- w_W*f; sumw <<- sumw*f }

d_cnt <- tabulate(cty[delta == 1], nbins = N)          # events per cluster (fixed)
sy    <- as.numeric(tapply(sd_df$time, factor(cty, levels = 1:N), sum))

for (iter in 1:num_iter) {
  ## --- data augmentation: rejected points of the thinned Poisson process ---
  rates <- lambda0 * Wv[cty] * sd_df$time
  q <- rpois(nt, pmax(rates, 1e-8))
  idx <- which(q > 0)
  Xaug <- NULL; m_cnt <- rep(0, N)
  if (length(idx) > 0) {
    rep_i <- rep(idx, q[idx]); G <- runif(length(rep_i), 0, sd_df$time[rep_i])
    XG <- cbind(G/max_time, M_hat[cty[rep_i]], X[rep_i, , drop = FALSE])
    keepG <- runif(nrow(XG)) > pnorm(forest$do_predict(XG))
    if (any(keepG)) {
      Xaug <- XG[keepG, , drop = FALSE]
      m_cnt <- tabulate(cty[rep_i[keepG]], nbins = N)
    }
  }
  ev_idx <- which(delta == 1)
  Xall <- rbind(Xaug, Xd[ev_idx, , drop = FALSE])
  Yall <- c(rep(0, if (is.null(Xaug)) 0 else nrow(Xaug)), rep(1, length(ev_idx)))
  bcur <- forest$do_predict(Xall)
  Z <- ifelse(Yall == 1, rtruncnorm(length(Yall), a = 0, b = Inf, mean = bcur, sd = 1),
              rtruncnorm(length(Yall), a = -Inf, b = 0, mean = bcur, sd = 1))
  forest$do_gibbs(Xall, Z, Xall, 1)
  
  ## --- lambda0 | . (Gibbs, Gamma(1,1) prior) ---
  lambda0 <- rgamma(1, 1 + sum(delta) + sum(m_cnt), 1 + sum(Wv[cty] * sd_df$time))
  ## --- sigma1^2 | . (Gibbs, IG(1,1) prior) ---
  qf <- as.numeric(t(Rv) %*% (D - rho1 * A) %*% Rv)
  s1 <- 1/rgamma(1, 1 + N/2, 1 + qf/2)
  ## --- rho1 | . (uniform random-walk MH) ---
  rp <- runif(1, max(rho_bounds[1], rho1 - 0.05), min(rho_bounds[2], rho1 + 0.05))
  if (log(runif(1)) < car_ld(Rv, rp, s1) - car_ld(Rv, rho1, s1)) { rho1 <- rp; acc_rho1 <- acc_rho1 + 1 }
  ## --- R = log W | . (component-wise Gaussian random-walk MH) ---
  for (c2 in 1:N) {
    Rp <- Rv; Rp[c2] <- rnorm(1, Rv[c2], 0.2)
    llp <- (d_cnt[c2] + m_cnt[c2]) * Rp[c2] - lambda0 * exp(Rp[c2]) * sy[c2]
    llc <- (d_cnt[c2] + m_cnt[c2]) * Rv[c2] - lambda0 * exp(Rv[c2]) * sy[c2]
    la <- (llp + car_ld(Rp, rho1, s1)) - (llc + car_ld(Rv, rho1, s1))
    if (is.finite(la) && log(runif(1)) < la) { Rv <- Rp; acc_R[c2] <- acc_R[c2] + 1 }
  }
  Wv <- exp(Rv)      # FIX 2: no deterministic recentring (see header)
  keep_lam[iter] <- lambda0; keep_s1[iter] <- s1; keep_rho1[iter] <- rho1; keep_W[iter, ] <- Wv
  
  ## --- posterior-mean survival accumulation (POST burn-in only) ---
  if (iter > burn_in && iter %% THIN == 0) {
    lw <- 0
    if (use_step3 == 1) {
      Mdraw <- M_post[sample.int(nrow(M_post), 1), ]     # draw M from p(M | D_0)
      Xq_m <- Xq_base; Xq_m[, 2] <- Mdraw[cty[q_rows]]
      Xq_h <- Xq_base; Xq_h[, 2] <- M_hat[cty[q_rows]]
      Xy_m <- Xy_base; Xy_m[, 2] <- Mdraw[cty]
      Xy_h <- Xy_base; Xy_h[, 2] <- M_hat[cty]
      Im  <- rowMeans(matrix(pnorm(forest$do_predict(Xq_m)), nt, NQ)) * sd_df$time
      Ih  <- rowMeans(matrix(pnorm(forest$do_predict(Xq_h)), nt, NQ)) * sd_df$time
      lPm <- log(pmax(pnorm(forest$do_predict(Xy_m)), 1e-300))
      lPh <- log(pmax(pnorm(forest$do_predict(Xy_h)), 1e-300))
      lw  <- sum(delta * (lPm - lPh)) - lambda0 * sum(Wv[cty] * (Im - Ih))
      lw_all <- c(lw_all, lw)
    }
    if (lw > maxlw) { if (is.finite(maxlw)) rescale_all(exp(maxlw - lw)); maxlw <- lw; w <- 1
    } else w <- exp(lw - maxlw)
    n_acc <- n_acc + 1; sumw <- sumw + w
    w_lam <- w_lam + w * lambda0; w_W <- w_W + w * Wv
    dt <- tg[1]
    clip01 <- function(m) { m[m < 0] <- 0; m[m > 1] <- 1; m }
    for (i in 1:N) {
      ph_i <- matrix(pnorm(forest$do_predict(amse_blocks[[i]])), nXP, nT0)
      S_amse_sum[i, , ] <- S_amse_sum[i, , ] +
        w * clip01(exp(-lambda0 * Wv[i] * cumtrap(ph_i, dt)))
    }
    ph_a <- matrix(pnorm(forest$do_predict(aes_block)), nPT0, 4)
    S_a  <- clip01(exp(-lambda0 * Wv[1] * cumtrap(t(ph_a), dtp)))
    S_aes_sum <- S_aes_sum + w * S_a
    ## store the raw draw and its log weight for the credible band
    n_keep <- n_keep + 1
    S_aes_draws[, , n_keep] <- S_a; lw_keep[n_keep] <- lw
    
    if (use_step3 == 1) {
      aesM <- aesM_base; aesM[, 2] <- Mdraw[1]
      ph_aM <- matrix(pnorm(forest$do_predict(aesM)), nPT0, 4)
      S_aesM_draws[, , n_keep] <- clip01(exp(-lambda0 * Wv[1] * cumtrap(t(ph_aM), dtp)))
    }
  }
  if (iter %% 1000 == 0) lg(sprintf("iter %5d  lam0=%.3f s1=%.2f rho1=%.2f accR=%.2f accRho=%.2f (%.1f min)",
                                    iter, lambda0, s1, rho1, mean(acc_R)/iter, acc_rho1/iter,
                                    as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

S_sb_amse <- S_amse_sum/sumw; S_sb_aes <- S_aes_sum/sumw
W_hat <- if (use_step3 == 1) w_W/sumw else colMeans(keep_W[(burn_in+1):num_iter, ])
ess   <- if (use_step3 == 1 && length(lw_all) > 0) {
  lwc <- lw_all - max(lw_all); ww <- exp(lwc); sum(ww)^2/sum(ww^2) } else n_acc

################################################################################
# 7. Metrics
################################################################################
amse_of <- function(Sest, k = 1:nT)
  mean(sapply(1:N, function(i) mean((S_true_amse[i,,k] - Sest[i,,k])^2)))
idx95 <- which(tg <= quantile(sd_df$time, 0.95))    # horizon used by the released code
amse_sbart <- amse_of(S_sb_amse)

# RSF / PH on the same grids
xfull <- cbind(xp, matrix(0.5, nXP, p - 3))
S_rsf_amse <- array(NA, c(N, nXP, nT)); S_ph_amse <- array(NA, c(N, nXP, nT))
for (i in 1:N) {
  if (!is.null(rsf)) S_rsf_amse[i,,] <- rsf_S(xfull, rep(Mplug[i], nXP), tg)
  if (!is.null(ph)) {
    lp <- as.numeric(xfull %*% bph[paste0("x",1:p)]) + bph["M"] * Mplug[i]
    S_ph_amse[i,,] <- exp(-W_ph[i] * outer(exp(lp), H0(tg)))
  }
}
amse_rsf <- if (!is.null(rsf)) amse_of(S_rsf_amse) else NA
amse_ph  <- if (!is.null(ph))  amse_of(S_ph_amse)  else NA

xaes <- t(sapply(aes_cases, function(a) c(a["x1"], 0.5, 0.5, rep(0.5, p - 3))))
S_rsf_aes <- if (!is.null(rsf)) rsf_S(xaes, sapply(aes_cases, function(a) a["p0"]), pt) else matrix(NA,4,nPT)
S_ph_aes  <- if (!is.null(ph)) {
  lp <- as.numeric(xaes %*% bph[paste0("x",1:p)]) + bph["M"] * sapply(aes_cases, function(a) a["p0"])
  exp(-W_ph[1] * outer(exp(lp), H0(pt))) } else matrix(NA, 4, nPT)


lg("intervals: SBART weighted credible band")
lwk   <- lw_keep[1:n_keep]
Wt    <- exp(lwk - max(lwk))               # Step-3 weights; all 1 when step3 = 0
ess_w <- sum(Wt)^2 / sum(Wt^2)             # what certifies the weighted band

wquant <- function(v, w, probs = c(0.025, 0.975)) {
  o <- order(v); cw <- cumsum(w[o])/sum(w)
  vapply(probs, function(pp) v[o][which(cw >= pp)[1]], numeric(1)) }
qq <- apply(S_aes_draws[, , 1:n_keep, drop = FALSE], c(1, 2), wquant, w = Wt)
sb_lo <- qq[1, , ]; sb_hi <- qq[2, , ]                    # M_1 fixed at the panel p0
sbM_lo <- sbM_hi <- NULL
if (use_step3 == 1) {                                     # M_1 ~ p(M | D_0) as well
  qqM <- apply(S_aesM_draws[, , 1:n_keep, drop = FALSE], c(1, 2), wquant, w = Wt)
  sbM_lo <- qqM[1, , ]; sbM_hi <- qqM[2, , ] }
lg(sprintf("  %d retained draws, importance ESS = %.1f (%.0f%%)",
           n_keep, ess_w, 100*ess_w/n_keep))

lg("intervals: PH-frailty bootstrap")
NB_PH <- 500
gridX <- as.data.frame(xaes); names(gridX) <- paste0("x", 1:p)
gridX$M <- as.numeric(sapply(aes_cases, function(a) a["p0"]))
Xg <- as.matrix(gridX[, paste0("x", 1:p)])
ph_draws <- array(NA_real_, c(4, nPT, NB_PH))
for (b in 1:NB_PH) {
  db <- ph_dat[sample.int(nt, nt, replace = TRUE), ]
  fb <- tryCatch(coxph(ph_fml, data = db), error = function(e) NULL)
  if (is.null(fb)) next
  bb <- coef(fb); bb[is.na(bb)] <- 0
  frb  <- fb$frail
  frbv <- if (!is.null(frb) && length(frb) == N) as.numeric(frb) else rep(0, N)
  lpb <- as.numeric(as.matrix(db[, paste0("x", 1:p)]) %*% bb[paste0("x", 1:p)]) +
    bb["M"] * db$M + frbv[db$county]
  ob <- order(db$time); tb <- db$time[ob]; dbo <- db$delta[ob]; rb <- exp(lpb)[ob]
  arb <- rev(cumsum(rev(rb))); etb <- unique(tb[dbo == 1])
  if (length(etb) < 2) next
  dkb <- as.numeric(table(factor(tb[dbo == 1], levels = etb)))
  H0b <- approx(etb, cumsum(dkb / arb[match(etb, tb)]), xout = pt,
                method = "constant", rule = 2, f = 0, ties = "ordered")$y
  lpg <- as.numeric(Xg %*% bb[paste0("x", 1:p)]) + bb["M"] * gridX$M
  ph_draws[, , b] <- exp(-exp(frbv[1]) * outer(exp(lpg), H0b))
  if (b %% 100 == 0) lg(sprintf("  PH bootstrap %d/%d", b, NB_PH))
}
ph_lo <- apply(ph_draws, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
ph_hi <- apply(ph_draws, c(1,2), quantile, probs = 0.975, na.rm = TRUE)

lg("intervals: RSF bootstrap")
NB_RSF <- 500; NT_B <- 250
rsf_draws <- array(NA_real_, c(4, nPT, NB_RSF))
for (b in 1:NB_RSF) {
  fb <- tryCatch(ranger(Surv(time, delta) ~ ., data = rsf_df[sample.int(nt, nt, replace = TRUE), ],
                        num.trees = NT_B, min.node.size = 15, seed = 1000*rep_id + b,
                        num.threads = 1), error = function(e) NULL)
  if (is.null(fb)) next
  pb <- tryCatch(predict(fb, data = gridX, num.threads = 1), error = function(e) NULL)
  if (is.null(pb)) next
  rsf_draws[, , b] <- t(apply(pb$survival, 1, function(v)
    approx(pb$unique.death.times, v, xout = pt, rule = 2)$y))
  if (b %% 100 == 0) lg(sprintf("  RSF bootstrap %d/%d", b, NB_RSF))
}
rsf_lo <- apply(rsf_draws, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
rsf_hi <- apply(rsf_draws, c(1,2), quantile, probs = 0.975, na.rm = TRUE)

own_wid <- c(sbart = mean(sb_hi - sb_lo, na.rm = TRUE),
             sbart_Mint = if (is.null(sbM_lo)) NA_real_ else mean(sbM_hi - sbM_lo, na.rm = TRUE),
             ph = mean(ph_hi - ph_lo, na.rm = TRUE),
             rsf = mean(rsf_hi - rsf_lo, na.rm = TRUE))
lg(sprintf("interval widths: SBART=%.4f (M-integrated %.4f) PH=%.4f RSF=%.4f",
           own_wid["sbart"], own_wid["sbart_Mint"], own_wid["ph"], own_wid["rsf"]))

saveRDS(list(rep_id = rep_id, seed = 20260814 + rep_id, num_tree = num_tree,
            
             own_width = own_wid, n_keep = n_keep, ess_weights = ess_w,
             sb_lo = sb_lo, sb_hi = sb_hi, sbM_lo = sbM_lo, sbM_hi = sbM_hi,
             ph_lo = ph_lo, ph_hi = ph_hi, rsf_lo = rsf_lo, rsf_hi = rsf_hi, scenario = scenario,
             version = if (use_step3) 3 else 2, ess_step3 = ess, n_acc = n_acc,
             logw_sd = if (length(lw_all)) sd(lw_all) else NA,
             n = nt, events = sum(delta), cens_frac = mean(delta == 0),
             amse = c(sbart = amse_sbart, rsf = amse_rsf, ph = amse_ph),
             amse_t95 = c(sbart = amse_of(S_sb_amse, idx95),
                          rsf = if (!is.null(rsf)) amse_of(S_rsf_amse, idx95) else NA,
                          ph  = if (!is.null(ph))  amse_of(S_ph_amse,  idx95) else NA),
             t95 = as.numeric(quantile(sd_df$time, 0.95)),
             cor_M_bayes = cor(M_true, M_hat), cor_M_naive = cor(M_true, Mplug),
             cor_W = cor(W_true, W_hat), W_true = W_true, W_hat = W_hat,
             M_true = M_true, M_hat = M_hat, M_plug = Mplug,
             rho1_post = mean(keep_rho1[(burn_in+1):num_iter]),
             s1_post = mean(keep_s1[(burn_in+1):num_iter]),
             lam0_post = mean(keep_lam[(burn_in+1):num_iter]),
             acc_R = mean(acc_R)/num_iter, acc_rho1 = acc_rho1/num_iter,
             plot_times = pt, S_true_aes = S_true_aes, S_sb_aes = S_sb_aes,
             S_rsf_aes = S_rsf_aes, S_ph_aes = S_ph_aes,
             trace = list(lam0 = keep_lam, s1 = keep_s1, rho1 = keep_rho1),
             minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
        file.path(outdir, sprintf("rep%02d.rds", rep_id)))

lg(sprintf("DONE  scen %s v%d  AMSE: SBART=%.5f RSF=%.5f PH=%.5f | ESS=%.1f | %.1f min",
           scenario, if (use_step3) 3 else 2, amse_sbart, amse_rsf, amse_ph, ess,
           as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
