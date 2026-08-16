#!/usr/bin/env Rscript
################################################################################
# Aggregate the Scenario A replicates: AMSE table + boxplot, AES figure with
# 95% pointwise bands, and the Table-S1 coverage row.
# Usage: Rscript scenA_aggregate.R <rundir>
################################################################################
rundir <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "scenA_v2_replication"
fs <- sort(list.files(rundir, pattern = "^rep[0-9]+\\.rds$", full.names = TRUE))
res <- lapply(fs, readRDS)
R <- length(res); cat("replicates:", R, "\n\n")

## ---------------- 1. AMSE ----------------------------------------------------
amse     <- t(sapply(res, function(x) x$amse))
amse95   <- t(sapply(res, function(x) x$amse_t95))
## published Scenario means/SDs and Table S1 coverage, selected by scenario
scen  <- if (!is.null(res[[1]]$scenario)) res[[1]]$scenario else "A"
vers  <- if (!is.null(res[[1]]$version))  res[[1]]$version  else 2
PM <- list(A = c(0.00381,0.01343,0.00798), B = c(0.00760,0.01704,0.02247),
           C = c(0.00534,0.01246,0.02581), D = c(0.00660,0.01210,0.07188))
PS <- list(A = c(0.00092,0.00269,0.00320), B = c(0.00208,0.00878,0.04003),
           C = c(0.00314,0.00235,0.05157), D = c(0.00080,0.00351,0.09037))
PC <- list(A = c(99.9,33.8,99.5), B = c(99.7,40.0,99.1),
           C = c(92.0,18.2,75.0), D = c(75.1,37.4,73.8))
paper_mu  <- setNames(PM[[scen]], c("sbart","rsf","ph"))
paper_sd  <- setNames(PS[[scen]], c("sbart","rsf","ph"))
paper_cvg <- setNames(PC[[scen]], c("sbart","rsf","ph"))
cat(sprintf("SCENARIO %s | version %d | %d replicates | K = %d trees\n\n", scen, vers, length(res),
            ifelse(is.null(res[[1]]$num_tree), 20, res[[1]]$num_tree)))

cat("=== AMSE, horizon (0, 10] as specified in the manuscript ===\n")
cat(sprintf("%-8s %10s %10s %10s | %10s %10s\n",
            "method", "mean", "SD", "median", "paper mean", "paper SD"))
for (m in c("sbart","rsf","ph"))
  cat(sprintf("%-8s %10.5f %10.5f %10.5f | %10.5f %10.5f\n", m,
              mean(amse[,m]), sd(amse[,m]), median(amse[,m]), paper_mu[m], paper_sd[m]))
cat("\n=== AMSE, horizon (0, t95] as used by the released code ===\n")
for (m in c("sbart","rsf","ph"))
  cat(sprintf("%-8s %10.5f %10.5f %10.5f\n", m,
              mean(amse95[,m]), sd(amse95[,m]), median(amse95[,m])))
cat(sprintf("\nSBART best in %d of %d replicates (both competitors)\n",
            sum(amse[,"sbart"] < amse[,"rsf"] & amse[,"sbart"] < amse[,"ph"]), R))

write.csv(data.frame(rep = sapply(res, `[[`, "rep_id"), amse,
                     amse_t95 = amse95), file.path(rundir, "amse_by_rep.csv"), row.names = FALSE)

## ---------------- 2. AES with 95% pointwise bands ----------------------------
pt   <- res[[1]]$plot_times
get  <- function(fld) simplify2array(lapply(res, `[[`, fld))      # 4 x nPT x R
Strue <- get("S_true_aes"); Ssb <- get("S_sb_aes")
Srsf  <- get("S_rsf_aes");  Sph <- get("S_ph_aes")
tq <- qt(0.975, R - 1)
clipm <- function(m, a, b) { m[m < a] <- a; m[m > b] <- b; m }
band <- function(Arr) {                                            # mean +/- t*SD
  m <- apply(Arr, c(1,2), mean, na.rm = TRUE); s <- apply(Arr, c(1,2), sd, na.rm = TRUE)
  list(m = m, lo = clipm(m - tq*s, 0, 1), hi = clipm(m + tq*s, 0, 1)) }
Bsb <- band(Ssb); Brsf <- band(Srsf); Bph <- band(Sph)
Btrue <- band(Strue)                    # the truth moves too: W_1 is redrawn each replicate
Tm  <- Btrue$m

labs <- c("x1=0.3,  M1=0.75", "x1=0.7,  M1=0.75", "x1=0.3,  M1=0.50", "x1=0.7,  M1=0.50")
draw <- function() {
  par(mfrow = c(2,2), mar = c(3.6,3.8,2.6,0.8), mgp = c(2.2,0.7,0))
  for (k in 1:4) {
    plot(pt, Tm[k,], type = "n", ylim = c(0,1), xlab = "t", ylab = "S(t)",
         main = labs[k], cex.main = 1)
    pb <- function(B, col) polygon(c(pt, rev(pt)), c(B$lo[k,], rev(B$hi[k,])),
                                   col = adjustcolor(col, 0.18), border = NA)
    polygon(c(pt, rev(pt)), c(Btrue$lo[k,], rev(Btrue$hi[k,])),
            col = adjustcolor("grey35", 0.16), border = NA)     # spread of the true curve
    pb(Bsb, "forestgreen"); pb(Brsf, "red"); pb(Bph, "blue")
    lines(pt, Tm[k,],     col = "black", lwd = 3.2)
    lines(pt, Bsb$m[k,],  col = "forestgreen", lwd = 2)
    lines(pt, Brsf$m[k,], col = "red",   lwd = 2, lty = 2)
    lines(pt, Bph$m[k,],  col = "blue",  lwd = 2, lty = 4)
    if (k == 1) legend("topright", c("mean true S (grey band = its spread)","spatial SBART","RSF","PH-frailty"),
                       col = c("black","forestgreen","red","blue"),
                       lty = c(1,1,2,4), lwd = 2, bty = "n", cex = 0.72)
  }
  mtext(sprintf("Scenario %s, version %d: AES at (x0, M1=p0) over %d replicates", scen, vers, R),
        outer = TRUE, line = 0.4, cex = 1.05)
}
pdf(file.path(rundir, paste0("AES_", scen, "_v", vers, ".pdf")), width = 8.5, height = 7); par(oma=c(0,0,2.4,0)); draw(); dev.off()
png(file.path(rundir, paste0("AES_", scen, "_v", vers, ".png")), width = 1700, height = 1400, res = 170); par(oma=c(0,0,2.4,0)); draw(); dev.off()

## ---------------- 3. AMSE boxplot -------------------------------------------
drawbox <- function() {
  par(mar = c(3.4,4.4,2.4,1))
  boxplot(list(`spatial SBART` = amse[,"sbart"], RSF = amse[,"rsf"], `PH-frailty` = amse[,"ph"]),
          col = c("grey95","grey80","grey55"), boxwex = 0.55, ylab = "AMSE",
          main = sprintf("Scenario %s, version %d: AMSE over %d replicates", scen, vers, R))
  points(1:3, paper_mu[c("sbart","rsf","ph")], pch = 4, cex = 1.6, lwd = 2, col = "red")
  legend("topleft", "x = AMSE mean reported in the paper", pch = 4, col = "red",
         bty = "n", cex = 0.85, pt.lwd = 2)
}
pdf(file.path(rundir, paste0("AMSE_", scen, "_v", vers, ".pdf")), width = 6, height = 5.2); drawbox(); dev.off()
png(file.path(rundir, paste0("AMSE_", scen, "_v", vers, ".png")), width = 1200, height = 1050, res = 170); drawbox(); dev.off()

## ---------------- 4. Coverage (Table S1, row A) ------------------------------
Tr_cov <- apply(simplify2array(lapply(res, `[[`, "S_true_cov")), 1:3, mean)
cov_of <- function(fld) {                    # (a) t-based pointwise band
  Est <- simplify2array(lapply(res, `[[`, fld))
  m <- apply(Est, 1:3, mean, na.rm = TRUE); s <- apply(Est, 1:3, sd, na.rm = TRUE)
  100 * mean(Tr_cov >= m - tq*s & Tr_cov <= m + tq*s, na.rm = TRUE) }
cov_emp <- function(fld) {                   # (b) empirical 2.5%/97.5% pointwise band
  Est <- simplify2array(lapply(res, `[[`, fld))
  lo <- apply(Est, 1:3, quantile, probs = 0.025, na.rm = TRUE)
  hi <- apply(Est, 1:3, quantile, probs = 0.975, na.rm = TRUE)
  100 * mean(Tr_cov >= lo & Tr_cov <= hi, na.rm = TRUE) }
cvg  <- c(sbart = cov_of("S_sb_cov"),  rsf = cov_of("S_rsf_cov"),  ph = cov_of("S_ph_cov"))
cvgE <- c(sbart = cov_emp("S_sb_cov"), rsf = cov_emp("S_rsf_cov"), ph = cov_emp("S_ph_cov"))
cat("\n=== Table S1 (Scenario A): % of 10,000 (t, x, M1) grid points covered ===\n")
cat(sprintf("%-8s %12s %12s %9s\n", "method", "t-band", "empirical", "paper"))
for (m in c("sbart","rsf","ph"))
  cat(sprintf("%-8s %11.1f%% %11.1f%% %8.1f%%\n", m, cvg[m], cvgE[m], paper_cvg[m]))

## ---------------- 5. Model-fit diagnostics ----------------------------------
cat("\n=== Estimation diagnostics (mean over replicates) ===\n")
d <- data.frame(
  cor_M_bayes = sapply(res, `[[`, "cor_M_bayes"), cor_M_naive = sapply(res, `[[`, "cor_M_naive"),
  cor_W = sapply(res, `[[`, "cor_W"), rho1 = sapply(res, `[[`, "rho1_post"),
  sigma1sq = sapply(res, `[[`, "s1_post"), lambda0 = sapply(res, `[[`, "lam0_post"),
  acc_R = sapply(res, `[[`, "acc_R"), acc_rho1 = sapply(res, `[[`, "acc_rho1"),
  cens = sapply(res, `[[`, "cens_frac"), minutes = sapply(res, `[[`, "minutes"))
print(round(colMeans(d), 4))
write.csv(cbind(rep = sapply(res, `[[`, "rep_id"), d), file.path(rundir, "diagnostics_by_rep.csv"), row.names = FALSE)
cat("\ntrue values: rho1 = 0.5, sigma1^2 = 1\n")
if (vers == 3) cat(sprintf("Step-3 weights: mean ESS %.1f of %.0f retained draws; mean sd(log omega) = %.1f\n",
    mean(sapply(res, function(x) x$ess_step3)), mean(sapply(res, function(x) x$n_acc)),
    mean(sapply(res, function(x) ifelse(is.null(x$logw_sd), NA, x$logw_sd)), na.rm = TRUE)))
cat("figures written to", rundir, "\n")
