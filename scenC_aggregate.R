rundir <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "scenA_out"
fs  <- sort(list.files(rundir, pattern = "^rep[0-9]+\\.rds$", full.names = TRUE))
res <- lapply(fs, readRDS)
R   <- length(res)

first <- function(f) { for (x in res) if (!is.null(x[[f]])) return(x[[f]]); NULL }
scen <- if (!is.null(first("scenario"))) first("scenario") else "A"

alg  <- if (!is.null(first("n_steps"))) first("n_steps") else if (!is.null(first("sbM_lo"))) 3L else 2L
ntree <- if (is.null(first("num_tree"))) NA_integer_ else first("num_tree")
cat(sprintf("SCENARIO %s | %d-step algorithm | %d replicate%s | K = %s trees\n",
            scen, alg, R, if (R == 1) "" else "s", ntree))
if (!is.null(first("settings"))) {
  s0 <- first("settings")
  cat(sprintf("settings: %d iterations, %d burn-in, THIN %d, NQ %d, fine grid %d,\n",
              s0["num_iter"], s0["burn_in"], s0["THIN"], s0["NQ"], s0["FINE_N"]))
  cat(sprintf("          bootstraps PH/RSF %d/%d, RSF trees fit/bootstrap %d/%d\n",
              s0["NB_PH"], s0["NB_RSF"], s0["RSF_TREES"], s0["NT_B"]))
}
cat("\n")

nm <- c(sbart = "spatial SBART", rsf = "RSF", ph = "PH-frailty")

amse   <- t(sapply(res, function(x) x$amse))
amse95 <- t(sapply(res, function(x) x$amse_t95))
if (R == 1) { amse <- matrix(amse, 1, dimnames = list(NULL, c("sbart","rsf","ph")))
              amse95 <- matrix(amse95, 1, dimnames = list(NULL, c("sbart","rsf","ph"))) }

q_rep <- function(v, p) if (R > 1) as.numeric(quantile(v, p, na.rm = TRUE)) else NA_real_

cat("=== AMSE of the estimated survival function (main text, Figure 2) ===\n")
cat("spread over the ", R, " replicate", if (R == 1) "" else "s",
    ", as empirical quantiles (no SD, no t)\n\n", sep = "")
cat(sprintf("%-15s %10s %10s %10s %10s\n", "method", "mean", "median", "2.5%", "97.5%"))
for (m in c("sbart","rsf","ph"))
  cat(sprintf("%-15s %10.5f %10.5f %10.5f %10.5f\n", nm[m], mean(amse[,m]),
              median(amse[,m]), q_rep(amse[,m], 0.025), q_rep(amse[,m], 0.975)))
cat("\n--- same, restricted to (0, t95] ---\n")
for (m in c("sbart","rsf","ph"))
  cat(sprintf("%-15s %10.5f %10.5f %10.5f %10.5f\n", nm[m], mean(amse95[,m]),
              median(amse95[,m]), q_rep(amse95[,m], 0.025), q_rep(amse95[,m], 0.975)))
cat(sprintf("\nspatial SBART has the lowest AMSE in %d of %d replicates\n",
            sum(amse[,"sbart"] < amse[,"rsf"] & amse[,"sbart"] < amse[,"ph"]), R))
write.csv(data.frame(rep = sapply(res, `[[`, "rep_id"), amse, amse_t95 = amse95),
          file.path(rundir, "amse_by_rep.csv"), row.names = FALSE)

pt   <- res[[1]]$plot_times
get  <- function(f) simplify2array(lapply(res, `[[`, f))
Strue <- get("S_true_aes"); Ssb <- get("S_sb_aes")
Srsf  <- get("S_rsf_aes");  Sph <- get("S_ph_aes")
if (R == 1) { d <- c(dim(res[[1]]$S_true_aes), 1)
              Strue <- array(Strue, d); Ssb <- array(Ssb, d)
              Srsf <- array(Srsf, d); Sph <- array(Sph, d) }

band <- function(Arr) {
  m  <- apply(Arr, c(1,2), mean, na.rm = TRUE)
  lo <- if (R > 1) apply(Arr, c(1,2), quantile, probs = 0.025, na.rm = TRUE) else m
  hi <- if (R > 1) apply(Arr, c(1,2), quantile, probs = 0.975, na.rm = TRUE) else m
  list(m = m, lo = lo, hi = hi) }
Bsb <- band(Ssb); Brsf <- band(Srsf); Bph <- band(Sph)
Tm  <- apply(Strue, c(1,2), mean, na.rm = TRUE)

labs <- c(expression(x[1]==0.3~","~~M[1]==0.75), expression(x[1]==0.7~","~~M[1]==0.75),
          expression(x[1]==0.3~","~~M[1]==0.50), expression(x[1]==0.7~","~~M[1]==0.50))
draw_aes <- function() {
  par(mfrow = c(2,2), mar = c(3.6,3.8,2.4,0.8), mgp = c(2.2,0.7,0))
  for (k in 1:4) {
    plot(pt, Tm[k,], type = "n", ylim = c(0,1), xlab = "t", ylab = "S(t)",
         main = labs[[k]], cex.main = 1)
    pb <- function(B, col) polygon(c(pt, rev(pt)), c(B$lo[k,], rev(B$hi[k,])),
                                   col = adjustcolor(col, 0.18), border = NA)
    if (R > 1) { pb(Bsb, "forestgreen"); pb(Brsf, "red"); pb(Bph, "blue") }
    lines(pt, Tm[k,],     col = "black",       lwd = 3)
    lines(pt, Bsb$m[k,],  col = "forestgreen", lwd = 2)
    lines(pt, Brsf$m[k,], col = "red",         lwd = 2, lty = 2)
    lines(pt, Bph$m[k,],  col = "blue",        lwd = 2, lty = 4)
    if (k == 1) legend("topright", c("true", "spatial SBART", "RSF", "PH-frailty"),
                       col = c("black","forestgreen","red","blue"),
                       lty = c(1,1,2,4), lwd = 2, bty = "n", cex = 0.8)
  }
}
fa <- file.path(rundir, sprintf("AES_%s", scen))
pdf(paste0(fa, ".pdf"), width = 8.5, height = 7); draw_aes(); invisible(dev.off())
png(paste0(fa, ".png"), width = 1700, height = 1400, res = 170); draw_aes(); invisible(dev.off())

draw_amse <- function() {
  par(mar = c(3.4,4.4,2.0,1))
  boxplot(list(`spatial SBART` = amse[,"sbart"], RSF = amse[,"rsf"], `PH-frailty` = amse[,"ph"]),
          col = c("grey95","grey80","grey55"), boxwex = 0.55, ylab = "AMSE",
          main = sprintf("Scenario %s (%d replicates)", scen, R))
}
fm <- file.path(rundir, sprintf("AMSE_%s", scen))
pdf(paste0(fm, ".pdf"), width = 6, height = 5); draw_amse(); invisible(dev.off())
png(paste0(fm, ".png"), width = 1200, height = 1000, res = 170); draw_amse(); invisible(dev.off())

if (!is.null(res[[1]]$own_width)) {
  wid <- c("sbart","sbart_Mint","rsf","ph")
  ow  <- t(sapply(res, function(x) x$own_width[wid]))
  if (R == 1) ow <- matrix(ow, 1, dimnames = list(NULL, wid))
  colnames(ow) <- wid
  st <- first("settings")
  nb <- if (is.null(st)) "" else sprintf(" (%d and %d resamples)", st["NB_PH"], st["NB_RSF"])
  cat("\n=== each method's own 95% pointwise interval on the AES grid ===\n")
  cat("SBART: posterior credible band (Step-3 weighted quantiles of the retained\n")
  cat(sprintf("draws).  PH, RSF: nonparametric bootstrap confidence bands%s.\n", nb))
  cat("No t distribution and no standard error is used; coverage is not computed.\n\n")
  cat(sprintf("%-26s %10s %10s %10s %10s\n", "mean width over the grid",
              "mean", "median", "2.5%", "97.5%"))
  for (m in wid) if (any(is.finite(ow[,m])))
    cat(sprintf("%-26s %10.4f %10.4f %10.4f %10.4f\n",
                if (m == "sbart") "SBART credible (M=p0)"
                else if (m == "sbart_Mint") "SBART credible (M drawn)"
                else if (m == "rsf") "RSF bootstrap" else "PH-frailty bootstrap",
                mean(ow[,m], na.rm = TRUE), median(ow[,m], na.rm = TRUE),
                q_rep(ow[,m], 0.025), q_rep(ow[,m], 0.975)))
  if (!is.null(res[[1]]$ess_weights)) {
    ess <- sapply(res, function(x) x$ess_weights); nk <- sapply(res, function(x) x$n_keep)
    cat(sprintf("\nimportance-sampling ESS behind the credible band: mean %.1f of %.0f retained draws (%.0f%%)\n",
                mean(ess), mean(nk), 100*mean(ess/nk)))
    cat("The ESS is what certifies the weighted band: the smaller it is, the fewer\n")
    cat("draws effectively support the 2.5% and 97.5% points.\n")
  }
  write.csv(data.frame(rep = sapply(res, `[[`, "rep_id"), width = ow,
                       ess = if (is.null(res[[1]]$ess_weights)) NA
                             else sapply(res, function(x) x$ess_weights),
                       n_keep = if (is.null(res[[1]]$n_keep)) NA
                                else sapply(res, function(x) x$n_keep)),
            file.path(rundir, "interval_width_by_rep.csv"), row.names = FALSE)

  grab <- function(f) { z <- lapply(res, `[[`, f)
                        if (any(sapply(z, is.null))) NULL else simplify2array(z) }
  saveRDS(list(plot_times = pt, aes_cases = list(c(x1=0.3,p0=0.75), c(x1=0.7,p0=0.75),
                                                 c(x1=0.3,p0=0.50), c(x1=0.7,p0=0.50)),
               rep_id = sapply(res, `[[`, "rep_id"),
               sb_lo = grab("sb_lo"), sb_hi = grab("sb_hi"),
               sbM_lo = grab("sbM_lo"), sbM_hi = grab("sbM_hi"),
               ph_lo = grab("ph_lo"), ph_hi = grab("ph_hi"),
               rsf_lo = grab("rsf_lo"), rsf_hi = grab("rsf_hi")),
          file.path(rundir, "intervals_by_rep.rds"))
  cat(sprintf("per-replicate bands written: %s\n",
              file.path(basename(rundir), "intervals_by_rep.rds")))
}
cat(sprintf("\nfigures written: %s.pdf|.png and %s.pdf|.png\n", basename(fa), basename(fm)))
