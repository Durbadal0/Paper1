# make_paper_outputs.R -- the simulation outputs of the paper, from the replicate files.
#   Rscript make_paper_outputs.R <results_dir> <out_dir>
# Writes to out_dir:
#   Figure2_AMSE_boxplots.pdf     Figure 2: per-replicate AMSE, one panel per scenario
#   FigureS1_AES_A.pdf ... S4     Figures S1-S4: true survival and AES curves; the band
#                                 is the 2.5%-97.5% range of the replicate AES curves
#   FigureS1_AES_A_ownband.pdf .. the same, with the band = average over replicates of
#                                 each method's own 95% pointwise interval
#   TableS1_coverage.csv / .tex   Table S1: share of the (t, x1, M1) grid points at which
#                                 each method's own 95% pointwise interval contains the
#                                 true survival function, averaged over replicates
#   amse_summary.csv              mean, SD, median, quartiles of AMSE per scenario/method
#   amse_by_replicate.csv         every replicate's AMSE
#   section4_numbers.txt          the numbers quoted in Section 4 and the run settings
args <- commandArgs(TRUE)
rdir <- if (length(args) >= 1) args[1] else "results"
odir <- if (length(args) >= 2) args[2] else "paper_outputs"
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

meth  <- c(sbart = "spatial SBART", rsf = "RSF", ph = "PH-frailty")
scens <- c("A", "B", "C", "D")
res <- list()
for (s in scens) {
  fs <- sort(list.files(file.path(rdir, paste0("scen", s)), "^rep[0-9]+\\.rds$", full.names = TRUE))
  if (length(fs)) res[[s]] <- lapply(fs, readRDS)
}
if (!length(res)) stop("no replicate files under ", rdir)
scens <- names(res)

# ---- AMSE
amse <- do.call(rbind, lapply(scens, function(s) {
  a <- t(sapply(res[[s]], function(x) x$amse[names(meth)]))
  data.frame(scenario = s, rep = sapply(res[[s]], `[[`, "rep_id"), a)
}))
write.csv(amse, file.path(odir, "amse_by_replicate.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(scens, function(s) do.call(rbind, lapply(names(meth), function(m) {
  v <- amse[amse$scenario == s, m]
  data.frame(scenario = s, method = meth[[m]], n_rep = sum(is.finite(v)),
             mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
             median = median(v, na.rm = TRUE), q25 = quantile(v, 0.25, na.rm = TRUE, names = FALSE),
             q75 = quantile(v, 0.75, na.rm = TRUE, names = FALSE),
             min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE))
}))))
write.csv(summ, file.path(odir, "amse_summary.csv"), row.names = FALSE)

# ---- Figure 2
draw_fig2 <- function() {
  par(mfrow = c(2, 2), mar = c(2.6, 4.4, 2.2, 0.8), mgp = c(2.6, 0.7, 0))
  for (s in scens) {
    a <- amse[amse$scenario == s, ]
    boxplot(list(`spatial SBART` = a$sbart, RSF = a$rsf, `PH-frailty` = a$ph),
            col = c("grey95", "grey75", "grey50"), boxwex = 0.55, ylab = "AMSE", las = 1,
            main = sprintf("Scenario %s", s), cex.axis = 0.85)
  }
}
pdf(file.path(odir, "Figure2_AMSE_boxplots.pdf"), width = 8, height = 6.5); draw_fig2(); invisible(dev.off())

# ---- Figures S1-S4
labs <- c(expression(x[1]==0.3~","~~M[1]==0.75), expression(x[1]==0.7~","~~M[1]==0.75),
          expression(x[1]==0.3~","~~M[1]==0.50), expression(x[1]==0.7~","~~M[1]==0.50))
stack <- function(s, f) simplify2array(lapply(res[[s]], `[[`, f))
draw_aes <- function(s, own) {
  pt <- res[[s]][[1]]$plot_times; R <- length(res[[s]])
  Tm <- apply(stack(s, "S_true_aes"), c(1, 2), mean)
  est <- list(sbart = stack(s, "S_sb_aes"), rsf = stack(s, "S_rsf_aes"), ph = stack(s, "S_ph_aes"))
  lo  <- list(sbart = "sb_lo", rsf = "rsf_lo", ph = "ph_lo"); hi <- list(sbart = "sb_hi", rsf = "rsf_hi", ph = "ph_hi")
  col <- c(sbart = "forestgreen", rsf = "red", ph = "blue"); lty <- c(sbart = 1, rsf = 2, ph = 4)
  par(mfrow = c(2, 2), mar = c(3.6, 3.8, 2.4, 0.8), mgp = c(2.2, 0.7, 0))
  for (k in 1:4) {
    plot(pt, Tm[k, ], type = "n", ylim = c(0, 1), xlab = "t", ylab = "S(t)", main = labs[[k]], cex.main = 1)
    for (m in names(meth)) {
      if (own) { L <- apply(stack(s, lo[[m]]), c(1, 2), mean, na.rm = TRUE)[k, ]
                 U <- apply(stack(s, hi[[m]]), c(1, 2), mean, na.rm = TRUE)[k, ] }
      else if (R > 1) { L <- apply(est[[m]], c(1, 2), quantile, 0.025, na.rm = TRUE)[k, ]
                        U <- apply(est[[m]], c(1, 2), quantile, 0.975, na.rm = TRUE)[k, ] }
      else next
      polygon(c(pt, rev(pt)), c(L, rev(U)), col = adjustcolor(col[[m]], 0.18), border = NA)
    }
    lines(pt, Tm[k, ], col = "black", lwd = 3)
    for (m in names(meth)) lines(pt, apply(est[[m]], c(1, 2), mean, na.rm = TRUE)[k, ],
                                 col = col[[m]], lwd = 2, lty = lty[[m]])
    if (k == 1) legend("topright", c("true", meth), col = c("black", col), lty = c(1, lty),
                       lwd = 2, bty = "n", cex = 0.8)
  }
}
for (i in seq_along(scens)) {
  s <- scens[i]; n <- match(s, c("A", "B", "C", "D"))
  pdf(file.path(odir, sprintf("FigureS%d_AES_%s.pdf", n, s)), width = 8.5, height = 7)
  draw_aes(s, own = FALSE); invisible(dev.off())
  pdf(file.path(odir, sprintf("FigureS%d_AES_%s_ownband.pdf", n, s)), width = 8.5, height = 7)
  draw_aes(s, own = TRUE); invisible(dev.off())
}

# ---- Table S1: coverage of each method's own 95% pointwise interval
covr <- do.call(rbind, lapply(scens, function(s) {
  v <- sapply(names(meth), function(m) {
    lo <- c(sbart = "sb_lo", rsf = "rsf_lo", ph = "ph_lo")[[m]]
    hi <- c(sbart = "sb_hi", rsf = "rsf_hi", ph = "ph_hi")[[m]]
    mean(sapply(res[[s]], function(x) mean(x[[lo]] <= x$S_true_aes & x$S_true_aes <= x[[hi]], na.rm = TRUE)))
  })
  data.frame(scenario = s, n_rep = length(res[[s]]), t(100 * v))
}))
write.csv(covr, file.path(odir, "TableS1_coverage.csv"), row.names = FALSE)
tex <- c("\\begin{tabular}{lccc}", "\\toprule", "Scenario & Spatial SBART & RSF & PH-frailty \\\\", "\\midrule",
         sprintf("%s & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\", covr$scenario, covr$sbart, covr$rsf, covr$ph),
         "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(odir, "TableS1_coverage.tex"))

# ---- numbers for Section 4
out <- c(sprintf("Simulation outputs from %s", normalizePath(rdir)), "")
for (s in scens) {
  x1 <- res[[s]][[1]]; st <- x1$settings
  ss <- summ[summ$scenario == s, ]
  a  <- amse[amse$scenario == s, ]
  out <- c(out, sprintf("Scenario %s: %d replicates; n = %d subjects; censoring %.1f%% (mean over replicates)",
                        s, length(res[[s]]), x1$n, 100 * mean(sapply(res[[s]], `[[`, "cens_frac"))),
           sprintf("  settings: %d iterations, %d burn-in, every %d-th kept (%d draws), K = %d trees, %d-step algorithm",
                   st[["num_iter"]], st[["burn_in"]], st[["THIN"]], x1$n_keep, st[["num_tree"]], x1$n_steps),
           sprintf("  bootstrap refits PH/RSF = %d/%d; run time per replicate: median %.0f min",
                   st[["NB_PH"]], st[["NB_RSF"]], median(sapply(res[[s]], `[[`, "minutes"))),
           sprintf("  AMSE mean (SD): %s", paste(sprintf("%s %.5f (%.5f)", ss$method, ss$mean, ss$sd), collapse = "; ")),
           sprintf("  AMSE median [IQR]: %s", paste(sprintf("%s %.5f [%.5f, %.5f]", ss$method, ss$median, ss$q25, ss$q75), collapse = "; ")),
           sprintf("  spatial SBART has the lowest AMSE in %d of %d replicates",
                   sum(a$sbart < a$rsf & a$sbart < a$ph, na.rm = TRUE), nrow(a)),
           sprintf("  coverage (Table S1): %s", paste(sprintf("%s %.1f%%", meth, unlist(covr[covr$scenario == s, names(meth)])), collapse = "; ")),
           "")
}
writeLines(out, file.path(odir, "section4_numbers.txt"))
cat(out, sep = "\n")
cat("outputs written to", normalizePath(odir), "\n")
