# sim_worker_cpp.R -- one replicate of the simulation study (Section 4), C++ version.
#
# The same script as sim_worker.R (the R version), except that the Step-1 loop and
# the Step-2/Step-3 loop are run by sbart_step1() and sbart_step2() in
# cpp/sbart_engine.cpp. Data generation, the RSF and PH-frailty fits, the AMSE and
# AES summaries and the bootstrap intervals are the same R code. The C++ loops draw
# their random numbers from R's generator in the same order as the R loops, so the
# two versions give the same results (see verify/).
#
# Usage: Rscript sim_worker_cpp.R rep outdir num_iter burn_in num_tree scenario step3 \
#                THIN NQ FINE_N NB_PH NB_RSF NT_B RSF_TREES NCORES
options(warn = -1)
# rep, outdir, iterations, burn-in, trees, scenario, step-3 switch
args <- commandArgs(trailingOnly = TRUE)
rep_id   <- as.integer(args[1])
outdir   <- if (length(args) >= 2) args[2] else "scenA_out"
num_iter <- if (length(args) >= 3) as.integer(args[3]) else 10000
burn_in  <- if (length(args) >= 4) as.integer(args[4]) else 3500
num_tree <- if (length(args) >= 5) as.integer(args[5]) else 200
scenario  <- if (length(args) >= 6) toupper(args[6]) else "A"
use_step3 <- if (length(args) >= 7) as.integer(args[7]) else 0

# optional overrides, args 8-14 (job.sh always passes them from config.sh)
argn <- function(k, d) if (length(args) >= k && nzchar(args[k])) as.integer(args[k]) else d
THIN      <- argn( 8, 3)
NQ        <- argn( 9, 50)
FINE_N    <- argn(10, 6001)
NB_PH     <- argn(11, 500)
NB_RSF    <- argn(12, 500)
NT_B      <- argn(13, 250)
RSF_TREES <- argn(14, 1000)
NCORES    <- argn(15, max(1L, as.integer(Sys.getenv("NCORES", "1"))))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({library(SoftBart); library(MASS); library(truncnorm)
                  library(survival); library(ranger)})
# compile (first time only) and load the C++ engine
local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  fa <- gsub("~+~", " ", sub("^--file=", "", fa), fixed = TRUE)   # Rscript writes spaces as ~+~
  code_dir <- if (length(fa)) dirname(normalizePath(fa[1])) else "."
  source(file.path(code_dir, "cpp", "load_engine.R"), local = TRUE)
  load_engine(code_dir)
})

logf <- file.path(outdir, sprintf("rep%02d.log", rep_id))
lg <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), " [rep ", rep_id, "] ", ...)
                      cat(m, "\n", file = logf, append = TRUE); cat(m, "\n"); flush.console() }
t_start <- Sys.time()

# the 15 panhandle counties and their real shared borders
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

# D is the small design: 10 counties, 1,000 subjects, t_max 8
if (scenario == "D") {
  A <- make_A10(); n_i <- c(80,100,120,100,80,120,80,100,120,100)
  n_0i <- c(8,12,15,20,10,18,9,25,14,11); t_max <- 8
} else {
  A <- make_A15(); n_i <- c(60,70,80,100,60,90,70,110,90,60,70,90,100,80,120)
  n_0i <- c(8,10,12,15,9,14,11,20,13,8,10,16,25,12,18); t_max <- 10
}
N <- nrow(A); D <- diag(rowSums(A))
p       <- 15
sig0    <- 1; sig1 <- 1; rho_true <- 0.5
cens_target <- 0.12

# admissible rho, from the normalised adjacency, capped at (-1,1)
dvec <- rowSums(A); Anorm <- diag(1/sqrt(dvec)) %*% A %*% diag(1/sqrt(dvec))
ev_n <- eigen((Anorm + t(Anorm))/2, only.values = TRUE)$values
rho_bounds <- c(max(-1, 1/min(ev_n)), min(1, 1/max(ev_n)))

# piecewise hazard; the pieces break at t = 2 and t = 6
g1 <- function(t,M,x1,x2,x3) 0.50*t*x1*x2 + 0.30*t*M*x1 + 0.15*M^2 + 0.40*x2*x3*M
g2 <- function(t,M,x1,x2,x3) 0.35*t*x1*x2*x3*M + 0.25*t*M^2 + 0.20*x1*M + 0.15*M^2
g3 <- function(t,M,x1,x2,x3) 0.30*log(t+1)*x1*x2*x3*M + 0.25*t*x2*x3 + 0.20*M^2

# D swaps the exponential link for a probit one
g4 <- function(t,M,x1,x2,x3) 0.5*t*x1*x2 + 0.3*t*M*x1 + 0.15*M^2 +
                             0.08*t^2*x1*x2*x3 + 0.3*t*x2*x3*M + 0.4*sqrt(x3)*M*t
haz <- function(t, W, M, x1, x2, x3) {
  if (scenario == "D") return(W * pnorm(g4(t,M,x1,x2,x3)))
  v <- ifelse(t < 2, g1(t,M,x1,x2,x3), ifelse(t < 6, g2(t,M,x1,x2,x3), g3(t,M,x1,x2,x3)))
  if (scenario == "B") { W^(1 + 0.2*t) * exp(v) } else { W * exp(v) } }

# true cumulative hazard, trapezoid on a fine grid
FINE <- seq(0, 30, length.out = FINE_N); dFINE <- FINE[2] - FINE[1]
cumH_fine <- function(W, M, x1, x2, x3) {
  h <- haz(FINE, W, M, x1, x2, x3)
  c(0, cumsum((h[-1] + h[-length(h)])/2 * dFINE)) }
S_true_fun <- function(tt, W, M, x1, x2, x3) {
  H <- cumH_fine(W, M, x1, x2, x3); exp(-approx(FINE, H, xout = tt, rule = 2)$y) }

# a fresh truth in every replicate
set.seed(20260814 + rep_id)
Qinv <- solve(D - rho_true * A + diag(1e-8, N))
R_true <- mvrnorm(1, rep(0, N), sig1^2 * Qinv); R_true <- R_true - mean(R_true)
W_true <- exp(R_true)

# C and D tie M to the frailty (spatial confounding)
M_true <- if (scenario %in% c("C","D")) {
  plogis(rnorm(N, 0.9 * R_true, sqrt(1 - 0.9^2)))
} else {
  plogis(as.numeric(mvrnorm(1, rep(0, N), sig0^2 * Qinv)))
}
m_0    <- rbinom(N, n_0i, M_true)
Mplug  <- m_0 / n_0i

nt <- sum(n_i); cty <- rep(1:N, times = n_i)
X  <- matrix(runif(nt * p), nt, p); colnames(X) <- paste0("x", 1:p)
Tobs <- numeric(nt)
for (k in 1:nt) {
  i <- cty[k]; E <- rexp(1)
  H <- cumH_fine(W_true[i], M_true[i], X[k,1], X[k,2], X[k,3])
  Tobs[k] <- if (max(H) < E) max(FINE) else approx(H, FINE, xout = E, ties = "ordered")$y
}

# censoring above the median, drawn independently of T
med   <- median(Tobs)
uexp  <- rexp(nt)
gap   <- function(r) mean(Tobs > med + uexp/r) - cens_target
r_cal <- tryCatch(uniroot(gap, c(1e-4, 1e4))$root, error = function(e) 1)
Cens  <- med + uexp/r_cal
delta <- as.integer(Tobs <= Cens)
Tobs  <- pmin(Tobs, Cens)
lg(sprintf("scenario %s | step3=%d | data: n=%d events=%d censored=%.1f%% median T=%.2f", scenario, use_step3,
           nt, sum(delta), 100*mean(delta == 0), med))
lg(sprintf("settings: K=%d THIN=%d NQ=%d FINE=%d bootstraps PH/RSF=%d/%d RSF trees fit/boot=%d/%d",
           num_tree, THIN, NQ, FINE_N, NB_PH, NB_RSF, RSF_TREES, NT_B))

sd_df <- data.frame(county = cty, time = Tobs, delta = delta, X)

# AMSE grid: 8x8x8 covariate points x 100 times
tg   <- seq(t_max/100, t_max, length.out = 100); nT <- length(tg)
gv   <- seq(0.05, 0.95, length.out = 8)
xp   <- as.matrix(expand.grid(x1 = gv, x2 = gv, x3 = gv)); nXP <- nrow(xp)
S_true_amse <- array(0, c(N, nXP, nT))
for (i in 1:N) for (r in 1:nXP)
  S_true_amse[i, r, ] <- S_true_fun(tg, W_true[i], M_true[i], xp[r,1], xp[r,2], xp[r,3])

# AES panels: county 1, four (x1, p0) combinations
t_aes <- 3
pt <- seq(t_aes/100, t_aes, length.out = 100); nPT <- length(pt); dtp <- pt[1]
aes_cases <- list(c(x1=0.3, p0=0.75), c(x1=0.7, p0=0.75),
                  c(x1=0.3, p0=0.50), c(x1=0.7, p0=0.50))
S_true_aes <- t(sapply(aes_cases, function(a)
  S_true_fun(pt, W_true[1], a["p0"], a["x1"], 0.5, 0.5)))

# competitor 1: RSF, with the naive plug-in for M
rsf_df <- cbind(sd_df[, c("time","delta")], as.data.frame(X), M = Mplug[cty])
rsf <- tryCatch(ranger(Surv(time, delta) ~ ., data = rsf_df, num.trees = RSF_TREES,
                       min.node.size = 15, seed = 42 + rep_id), error = function(e) NULL)
rsf_S <- function(xmat, Mv, tt) {
  nd <- as.data.frame(cbind(xmat, M = Mv)); names(nd) <- c(paste0("x",1:p), "M")
  pr <- predict(rsf, data = nd)
  t(apply(pr$survival, 1, function(s) approx(pr$unique.death.times, s, xout = tt, rule = 2)$y)) }

# competitor 2: PH frailty; sparse = FALSE keeps M identified
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

  # Breslow baseline, frailties kept in the risk set
  lp_dat <- as.numeric(as.matrix(ph_dat[, paste0("x", 1:p)]) %*% bph[paste0("x", 1:p)]) +
            bph["M"] * ph_dat$M + frailv[ph_dat$county]
  o   <- order(ph_dat$time)
  tso <- ph_dat$time[o]; dso <- ph_dat$delta[o]; rso <- exp(lp_dat)[o]
  atrisk <- rev(cumsum(rev(rso)))
  et  <- unique(tso[dso == 1])
  dk  <- as.numeric(table(factor(tso[dso == 1], levels = et)))
  den <- atrisk[match(et, tso)]
  H0_step <- cumsum(dk / den)
  H0 <- function(tt) approx(et, H0_step, xout = tt, method = "constant",
                            rule = 2, f = 0, ties = "ordered")$y
  lg(sprintf("PH fitted: beta_M = %.3f | %d event times | H0(max t) = %.3f",
             bph["M"], length(et), H0(max(sd_df$time))))
}

# ---- Step 1: CAR model for M from the survey counts
car_ld <- function(x, rho, s2) {
  Q <- (D - rho * A); ev <- eigen(Q, only.values = TRUE, symmetric = TRUE)$values
  0.5 * (sum(log(pmax(ev, 1e-12))) - N * log(s2)) - 0.5 * as.numeric(t(x) %*% Q %*% x) / s2 }
n1 <- 12000; b1 <- 6000
logitM <- qlogis((m_0 + 0.5)/(n_0i + 1)); s0 <- 1; rho0 <- 0
Msave <- matrix(NA, n1, N)
# the Step-1 loop of sim_worker.R, run in C++ (sbart_step1)
logdet_fun <- function(rho) {
  Q <- (D - rho * A); ev <- eigen(Q, only.values = TRUE, symmetric = TRUE)$values
  sum(log(pmax(ev, 1e-12))) }
# the quadratic form t(x) Q x, computed by R (it sets the sigma^2 draws)
qf_fun <- function(x, rho) as.numeric(t(x) %*% (D - rho * A) %*% x)
st1 <- sbart_step1(list(A = A, D = D, m_0 = as.integer(m_0), n_0i = as.numeric(n_0i),
                        logitM = logitM, s0 = s0, rho0 = rho0, rho_bounds = rho_bounds,
                        n1 = n1), logdet_fun, qf_fun)
Msave <- st1$Msave; logitM <- st1$logitM; s0 <- st1$s0; rho0 <- st1$rho0
M_post <- Msave[(b1+1):n1, ]
M_hat  <- colMeans(M_post)
lg(sprintf("step1 done: cor(M_true, M_hat) = %.3f  (naive plug-in cor = %.3f)",
           cor(M_true, M_hat), cor(M_true, Mplug)))

# ---- Step 2: SBART survival model
max_time <- max(sd_df$time)
Xd <- cbind(sd_df$time/max_time, M_hat[cty], X)
colnames(Xd) <- c("time", "M", paste0("x", 1:p))
Yinit <- ifelse(delta == 1, 0.5, -0.5) + rnorm(nt, 0, 0.1)

# sigma = 1 and sigma_mu = 3/(2 sqrt K), as the probit augmentation needs
hyp <- Hypers(Xd, Yinit, num_tree = num_tree, sigma_hat = 1, k = 1/3)
# the forest is built in C++ from these hypers and options (as MakeForest does)
opts_sb <- Opts(num_burn = 0, num_save = 1, update_sigma = FALSE)

lambda0 <- 0.1; Rv <- rep(0, N); Wv <- exp(Rv); s1 <- 1; rho1 <- 0
keep_lam <- keep_s1 <- keep_rho1 <- numeric(num_iter); keep_W <- matrix(NA, num_iter, N)
acc_R <- rep(0, N); acc_rho1 <- 0

xfill <- matrix(0.5, 1, p - 3)
tg0 <- c(0, tg); nT0 <- length(tg0)
cumtrap <- function(P, h) {
  m <- (P[, -1, drop = FALSE] + P[, -ncol(P), drop = FALSE]) / 2
  t(apply(m, 1, cumsum)) * h }
# the prediction blocks never change, so build them once
amse_blocks <- lapply(1:N, function(i) {
  bx <- cbind(xp, xfill[rep(1, nXP), , drop = FALSE])
  cbind(rep(tg0, each = nXP)/max_time, M_hat[i], bx[rep(1:nXP, times = nT0), ]) })
pt0 <- c(0, pt); nPT0 <- length(pt0)
aes_block <- do.call(rbind, lapply(aes_cases, function(a)
  cbind(pt0/max_time, a["p0"], a["x1"], 0.5, 0.5, xfill[rep(1, nPT0), , drop = FALSE])))

aesM_base <- aes_block

n_keep_max <- sum(seq_len(num_iter) > burn_in & seq_len(num_iter) %% THIN == 0)
S_amse_sum <- array(0, c(N, nXP, nT)); S_aes_sum <- matrix(0, 4, nPT); n_acc <- 0

# keep every retained AES draw; the credible band is built from these
S_aes_draws  <- array(NA_real_, c(4, nPT, n_keep_max))
S_aesM_draws <- if (use_step3 == 1) array(NA_real_, c(4, nPT, n_keep_max)) else NULL
lw_keep <- numeric(n_keep_max); n_keep <- 0

sumw <- 0; maxlw <- -Inf; lw_all <- c(); w_lam <- 0; w_W <- rep(0, N)
w_rho1 <- 0; w_s1 <- 0
# quadrature points for the Step-3 integral
qw <- (seq_len(NQ) - 0.5)/NQ
q_rows  <- rep(1:nt, times = NQ)
Xq_base <- cbind(as.numeric(outer(sd_df$time, qw))/max_time, 0, X[q_rows, , drop = FALSE])
Xy_base <- cbind(sd_df$time/max_time, 0, X)
# rescale the running sums whenever a bigger weight shows up
rescale_all <- function(f) { S_amse_sum <<- S_amse_sum*f; S_aes_sum <<- S_aes_sum*f
  w_lam <<- w_lam*f; w_W <<- w_W*f; w_rho1 <<- w_rho1*f; w_s1 <<- w_s1*f
  sumw <<- sumw*f }

d_cnt <- tabulate(cty[delta == 1], nbins = N)
sy    <- as.numeric(tapply(sd_df$time, factor(cty, levels = 1:N), sum))

# the Step-2/Step-3 loop of sim_worker.R, run in C++ (sbart_step2)
log_iter <- function(iter, lambda0, s1, rho1, acc_R, acc_rho1)
  lg(sprintf("iter %5d  lam0=%.3f s1=%.2f rho1=%.2f accR=%.2f accRho=%.2f (%.1f min)",
      iter, lambda0, s1, rho1, mean(acc_R)/iter, acc_rho1/iter,
      as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
st2 <- sbart_step2(list(X = X, Xd = Xd, A = A, D = D, time = sd_df$time,
                        delta = as.integer(delta), cty = as.integer(cty), M_hat = M_hat,
                        rho_bounds = rho_bounds, d_cnt = as.integer(d_cnt), sy = sy,
                        max_time = max_time, num_iter = num_iter, burn_in = burn_in,
                        THIN = THIN, use_step3 = use_step3, n_keep_max = n_keep_max,
                        nthreads = NCORES, amse_blocks = amse_blocks, aes_block = aes_block,
                        aesM_base = aesM_base, nXP = nXP, nT0 = nT0, nPT0 = nPT0,
                        dt = tg[1], dtp = dtp, M_post = M_post, NQ = NQ,
                        Xq_base = Xq_base, Xy_base = Xy_base,
                        cty_q = as.integer(cty[q_rows]), lambda0 = lambda0, s1 = s1,
                        rho1 = rho1, Rv = Rv),
                   hyp, opts_sb, logdet_fun, qf_fun, log_iter)
for (nm in c("keep_lam", "keep_s1", "keep_rho1", "keep_W", "acc_R", "acc_rho1",
             "S_amse_sum", "S_aes_sum", "sumw", "maxlw", "w_lam", "w_W", "w_rho1",
             "w_s1", "n_acc", "S_aes_draws", "lw_keep", "n_keep", "lambda0", "Rv",
             "Wv", "s1", "rho1"))
  assign(nm, st2[[nm]])
lw_all <- if (is.null(st2$lw_all)) c() else st2$lw_all
S_aesM_draws <- st2$S_aesM_draws
rm(st2)

S_sb_amse <- S_amse_sum/sumw; S_sb_aes <- S_aes_sum/sumw
W_hat <- if (use_step3 == 1) w_W/sumw else colMeans(keep_W[(burn_in+1):num_iter, ])
ess   <- if (use_step3 == 1 && length(lw_all) > 0) {
           lwc <- lw_all - max(lw_all); ww <- exp(lwc); sum(ww)^2/sum(ww^2) } else n_acc

amse_of <- function(Sest, k = 1:nT)
  mean(sapply(1:N, function(i) mean((S_true_amse[i,,k] - Sest[i,,k])^2)))
idx95 <- which(tg <= quantile(sd_df$time, 0.95))
amse_sbart <- amse_of(S_sb_amse)

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

# ---- each method's own 95% interval; no coverage is computed
lg("intervals: SBART weighted credible band")
lwk   <- lw_keep[1:n_keep]
Wt    <- exp(lwk - max(lwk))
ess_w <- sum(Wt)^2 / sum(Wt^2)

# weighted quantiles, so band and point estimate use the same draws
wquant <- function(v, w, probs = c(0.025, 0.975)) {
  o <- order(v); cw <- cumsum(w[o])/sum(w)
  vapply(probs, function(pp) v[o][which(cw >= pp)[1]], numeric(1)) }
qq <- apply(S_aes_draws[, , 1:n_keep, drop = FALSE], c(1, 2), wquant, w = Wt)
sb_lo <- qq[1, , ]; sb_hi <- qq[2, , ]
sbM_lo <- sbM_hi <- NULL
if (use_step3 == 1) {
  qqM <- apply(S_aesM_draws[, , 1:n_keep, drop = FALSE], c(1, 2), wquant, w = Wt)
  sbM_lo <- qqM[1, , ]; sbM_hi <- qqM[2, , ] }
lg(sprintf("  %d retained draws, importance ESS = %.1f (%.0f%%)",
           n_keep, ess_w, 100*ess_w/n_keep))


# run f(1), ..., f(n) on NCORES cores; f uses no random numbers, so the order of
# evaluation does not matter. A result lost by a crashed child is recomputed here.
par_map <- function(n, f) {
  if (NCORES <= 1 || n <= 1) return(lapply(seq_len(n), f))
  out <- parallel::mclapply(seq_len(n), f, mc.cores = NCORES, mc.preschedule = TRUE)
  bad <- which(!vapply(out, function(z) is.list(z) && isTRUE(z$done), logical(1)))
  if (length(bad)) { lg(sprintf("  %d refit(s) redone on one core", length(bad)))
                     out[bad] <- lapply(bad, f) }
  out }

# the main RSF forest (about 1.4 GB) and the SBART forest are not used any more;
# freeing them keeps the forked bootstrap workers small (no random numbers involved)
if (!is.null(rsf)) rsf <- "freed"   # a marker: the saved results test !is.null(rsf)
if (exists("forest")) rm(forest)
invisible(gc())

lg("intervals: PH-frailty bootstrap")
gridX <- as.data.frame(xaes); names(gridX) <- paste0("x", 1:p)
gridX$M <- as.numeric(sapply(aes_cases, function(a) a["p0"]))
Xg <- as.matrix(gridX[, paste0("x", 1:p)])
ph_draws <- array(NA_real_, c(4, nPT, NB_PH))
# every resample refits the whole frailty model, baseline included.
# The resamples are drawn here in the one-core order; the refits run on NCORES cores.
ph_idx <- lapply(1:NB_PH, function(b) sample.int(nt, nt, replace = TRUE))
ph_one <- function(b) {
  db <- ph_dat[ph_idx[[b]], ]
  fb <- tryCatch(coxph(ph_fml, data = db), error = function(e) NULL)
  if (is.null(fb)) return(list(done = TRUE, S = NULL))
  bb <- coef(fb); bb[is.na(bb)] <- 0
  frb  <- fb$frail
  frbv <- if (!is.null(frb) && length(frb) == N) as.numeric(frb) else rep(0, N)
  lpb <- as.numeric(as.matrix(db[, paste0("x", 1:p)]) %*% bb[paste0("x", 1:p)]) +
         bb["M"] * db$M + frbv[db$county]
  ob <- order(db$time); tb <- db$time[ob]; dbo <- db$delta[ob]; rb <- exp(lpb)[ob]
  arb <- rev(cumsum(rev(rb))); etb <- unique(tb[dbo == 1])
  if (length(etb) < 2) return(list(done = TRUE, S = NULL))
  dkb <- as.numeric(table(factor(tb[dbo == 1], levels = etb)))
  H0b <- approx(etb, cumsum(dkb / arb[match(etb, tb)]), xout = pt,
                method = "constant", rule = 2, f = 0, ties = "ordered")$y
  lpg <- as.numeric(Xg %*% bb[paste0("x", 1:p)]) + bb["M"] * gridX$M
  list(done = TRUE, S = exp(-exp(frbv[1]) * outer(exp(lpg), H0b)))
}
ph_res <- par_map(NB_PH, ph_one)
for (b in 1:NB_PH) if (!is.null(ph_res[[b]]$S)) ph_draws[, , b] <- ph_res[[b]]$S
lg(sprintf("  PH bootstrap %d/%d (%d cores)", NB_PH, NB_PH, NCORES))
ph_lo <- apply(ph_draws, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
ph_hi <- apply(ph_draws, c(1,2), quantile, probs = 0.975, na.rm = TRUE)

lg("intervals: RSF bootstrap")
rsf_draws <- array(NA_real_, c(4, nPT, NB_RSF))
# The one-core loop draws, for each b, the resample and then the runif() that
# predict() makes when it is given no seed. Both are drawn here in that order and
# passed to the refits, which run on NCORES cores.
rng_before_rsf <- .Random.seed
rsf_idx <- vector("list", NB_RSF); rsf_pseed <- numeric(NB_RSF)
for (b in 1:NB_RSF) { rsf_idx[[b]] <- sample.int(nt, nt, replace = TRUE)
                      rsf_pseed[b] <- runif(1, 0, .Machine$integer.max) }
rsf_one <- function(b) {
  fb <- tryCatch(ranger(Surv(time, delta) ~ ., data = rsf_df[rsf_idx[[b]], ],
                        num.trees = NT_B, min.node.size = 15, seed = 1000*rep_id + b,
                        num.threads = 1), error = function(e) NULL)
  if (is.null(fb)) return(list(done = TRUE, ok = FALSE, S = NULL))
  pb <- tryCatch(predict(fb, data = gridX, num.threads = 1, seed = rsf_pseed[b]),
                 error = function(e) NULL)
  if (is.null(pb)) return(list(done = TRUE, ok = FALSE, S = NULL))
  list(done = TRUE, ok = TRUE, S = t(apply(pb$survival, 1, function(v)
                        approx(pb$unique.death.times, v, xout = pt, rule = 2)$y)))
}
rsf_res <- par_map(NB_RSF, rsf_one)
if (all(vapply(rsf_res, function(z) isTRUE(z$ok), logical(1)))) {
  for (b in 1:NB_RSF) rsf_draws[, , b] <- rsf_res[[b]]$S
  lg(sprintf("  RSF bootstrap %d/%d (%d cores)", NB_RSF, NB_RSF, NCORES))
} else {
  # a failed refit skips its runif() in the one-core loop, so rerun that loop as written
  lg("  a refit failed: RSF bootstrap rerun on one core, in the original order")
  assign(".Random.seed", rng_before_rsf, envir = globalenv())
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
}
rsf_lo <- apply(rsf_draws, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
rsf_hi <- apply(rsf_draws, c(1,2), quantile, probs = 0.975, na.rm = TRUE)

own_wid <- c(sbart = mean(sb_hi - sb_lo, na.rm = TRUE),
             sbart_Mint = if (is.null(sbM_lo)) NA_real_ else mean(sbM_hi - sbM_lo, na.rm = TRUE),
             ph = mean(ph_hi - ph_lo, na.rm = TRUE),
             rsf = mean(rsf_hi - rsf_lo, na.rm = TRUE))
lg(sprintf("interval widths: SBART=%.4f (M-integrated %.4f) PH=%.4f RSF=%.4f",
           own_wid["sbart"], own_wid["sbart_Mint"], own_wid["ph"], own_wid["rsf"]))

# everything this replicate produced
# written to a .part file and renamed when complete, so that a job killed while
# saving leaves no half-written file that would count as finished
rds_f <- file.path(outdir, sprintf("rep%02d.rds", rep_id))
saveRDS(list(rep_id = rep_id, seed = 20260814 + rep_id, num_tree = num_tree,

  settings = c(num_iter = num_iter, burn_in = burn_in, num_tree = num_tree,
               THIN = THIN, NQ = NQ, FINE_N = FINE_N, NB_PH = NB_PH,
               NB_RSF = NB_RSF, NT_B = NT_B, RSF_TREES = RSF_TREES),

  own_width = own_wid, n_keep = n_keep, ess_weights = ess_w, logw = lwk,
  sb_lo = sb_lo, sb_hi = sb_hi, sbM_lo = sbM_lo, sbM_hi = sbM_hi,
  ph_lo = ph_lo, ph_hi = ph_hi, rsf_lo = rsf_lo, rsf_hi = rsf_hi, scenario = scenario,
  n_steps = if (use_step3) 3L else 2L, ess_step3 = ess, n_acc = n_acc,
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
  rho1_post = if (use_step3 == 1) w_rho1/sumw else mean(keep_rho1[(burn_in+1):num_iter]),
  s1_post = if (use_step3 == 1) w_s1/sumw else mean(keep_s1[(burn_in+1):num_iter]),
  lam0_post = if (use_step3 == 1) w_lam/sumw else mean(keep_lam[(burn_in+1):num_iter]),
  acc_R = mean(acc_R)/num_iter, acc_rho1 = acc_rho1/num_iter,
  plot_times = pt, S_true_aes = S_true_aes, S_sb_aes = S_sb_aes,
  S_rsf_aes = S_rsf_aes, S_ph_aes = S_ph_aes,
  trace = list(lam0 = keep_lam, s1 = keep_s1, rho1 = keep_rho1),
  minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
  paste0(rds_f, ".part"))
invisible(file.rename(paste0(rds_f, ".part"), rds_f))

lg(sprintf("DONE  scen %s (%d-step)  AMSE: SBART=%.5f RSF=%.5f PH=%.5f | ESS=%.1f | %.1f min",
           scenario, if (use_step3) 3L else 2L, amse_sbart, amse_rsf, amse_ph, ess,
           as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
