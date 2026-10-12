#!/usr/bin/env Rscript
################################################################################
# Tables and figures of the Data Analysis section (main text and supplement).
#
# Reads the outputs of fcr_data_analysis.R, one folder per stratum:
#   <outdir>/<RACE>_<STAGE>/fcr_<RACE>_<STAGE>_details.rds
#   (optional) <outdir>/AA_D_chain2/fcr_AA_D_details.rds  -- second chain, Figure S17
# and the data files in the working directory (not distributed):
#   Surv_data2.csv, mammo_county_16A.xlsx, mammo_county_16W.xlsx
# Strata whose output is missing are skipped (their columns/curves are left out).
#
# Writes to <figdir>/:
#   logrank.txt                          stratified log-rank test (Section 5)
#   table1_vip.csv / .tex                Table 1
#   figure3_boxplots.pdf                 Figure 3
#   figure4_survival.pdf                 Figure 4
#   table2_dms10.csv / .tex              Table 2 (and the DMS values quoted in the text)
#   figure5_dms10_smu.pdf                Figure 5
#   figS5_km.pdf                         Figure S5
#   figS6_coxsnell.pdf                   Figure S6
#   figS7_deviance_AA_R.pdf              Figure S7
#   figS8_deviance_WA_R.pdf              Figure S8 (shares outside +-2 in residual_shares.csv)
#   tableS2_dms5.csv / .tex              Table S2
#   figS9_map_dms10_smu.pdf              Figure S9  (needs the 'maps' package)
#   figS10_map_dms10_td.pdf              Figure S10 (needs the 'maps' package)
#   figS11_dms10_td.pdf                  Figure S11
#   figS12_dms10_td_ci.pdf               Figure S12
#   figS13_var_logitM.pdf                Figure S13
#   figS14_Mhat_comparison.pdf           Figure S14
#   figS15_Mhat_ci.pdf                   Figure S15
#   figS16_W_ci.pdf                      Figure S16
#   figS17_traces.pdf                    Figure S17
#   county_dms10_AA_R.csv                per-county values behind Figures 5 and S9-S12
#
# Usage: Rscript make_figures_tables.R [outdir = out] [figdir = figures]
################################################################################
suppressMessages({ library(survival); library(readxl) })
args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "out"
figdir <- if (length(args) >= 2) args[2] else "figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fp <- function(f) file.path(figdir, f)

## Coding of the delay variables. fcr_data_analysis.R hard-codes TX_Delay == 3 and
## BX_Delay == 3 as the indicators of a SHORT delay; keep SHORT_CODE equal to that.
SHORT_CODE <- 3

## the 67 Florida counties in the order of the county codes of the data
## (code = (FIPS + 1)/2, with Miami-Dade at its former FIPS code 025, i.e. county 13)
COUNTY <- c("Alachua", "Baker", "Bay", "Bradford", "Brevard", "Broward", "Calhoun",
  "Charlotte", "Citrus", "Clay", "Collier", "Columbia", "Miami-Dade", "DeSoto", "Dixie",
  "Duval", "Escambia", "Flagler", "Franklin", "Gadsden", "Gilchrist", "Glades", "Gulf",
  "Hamilton", "Hardee", "Hendry", "Hernando", "Highlands", "Hillsborough", "Holmes",
  "Indian River", "Jackson", "Jefferson", "Lafayette", "Lake", "Lee", "Leon", "Levy",
  "Liberty", "Madison", "Manatee", "Marion", "Martin", "Monroe", "Nassau", "Okaloosa",
  "Okeechobee", "Orange", "Osceola", "Palm Beach", "Pasco", "Pinellas", "Polk", "Putnam",
  "St. Johns", "St. Lucie", "Santa Rosa", "Sarasota", "Seminole", "Sumter", "Suwannee",
  "Taylor", "Union", "Volusia", "Wakulla", "Walton", "Washington")
N <- length(COUNTY)
## counties highlighted in Figures 5 and S11 (symbols as in the captions)
HIGHLIGHT <- c(Baker = 2, Escambia = 17, Broward = 6, Citrus = 9)
HL_PCH    <- c(Baker = 1, Escambia = 0, Broward = 8, Citrus = 2)
## counties whose log-frailty traces are shown in Figure S17 (Broward, Citrus, Duval)
TRACE_COUNTIES <- c(6, 9, 16)

################################################################################
# Outputs of fcr_data_analysis.R
################################################################################
STRATA <- c("AA_L", "WA_L", "AA_R", "WA_R", "AA_D", "WA_D")       # Table 1 order
load_run <- function(key, sub = key) {
  f <- file.path(outdir, sub, sprintf("fcr_%s_details.rds", key))
  if (file.exists(f)) readRDS(f) else NULL
}
runs <- setNames(lapply(STRATA, load_run), STRATA)
have <- STRATA[!vapply(runs, is.null, logical(1))]
cat("strata found:", if (length(have)) paste(have, collapse = ", ") else "none", "\n")
if (!length(have)) stop("no outputs of fcr_data_analysis.R found under ", outdir)
for (k in have) { st <- runs[[k]]$settings
  if (st[["use_step3"]] != 1) stop(k, " was run without Step 3 (use_step3 = 0); all tables and figures need the 3-step method")
  if (st[["spatial"]] != 1 || st[["t_hor"]] != 10) stop(k, ": the tables and figures need spatial = 1 and t_horizon = 10") }
missing_note <- function(keys) {
  m <- setdiff(keys, have)
  if (length(m)) cat("  (skipped, no output yet:", paste(m, collapse = ", "), ")\n")
}
f3 <- function(v) sprintf("%.3f", round(v, 3) + 0)        # 3 decimals, no "-0.000"

################################################################################
# Data, processed exactly as in fcr_data_analysis.R
################################################################################
surv_raw <- read.csv("Surv_data2.csv")
names(surv_raw) <- c("id", "time_days", "county", "TX_Delay", "death",
                     "BX_Delay", "Age", "HR_p", "Tgrade", "Race", "SMUprob", "Stage")
surv_raw$yr    <- pmax(surv_raw$time_days/365.25, 1/365.25)
surv_raw$delta <- as.integer(surv_raw$death == 1)
surv_raw$race  <- factor(ifelse(surv_raw$Race == 2, "AA", "WA"), levels = c("AA", "WA"))
surv_raw$stage <- factor(c("Local", "Regional", "Distant")[surv_raw$Stage],
                         levels = c("Local", "Regional", "Distant"))

stratum_data <- function(race, stage) {
  st <- surv_raw[surv_raw$Race == (if (race == "AA") 2 else 1) &
                 surv_raw$Stage == switch(stage, L = 1, R = 2, D = 3), ]
  d <- data.frame(county = st$county, time = st$yr, delta = st$delta, Age = st$Age,
                  HR = as.integer(st$HR_p == 1), Grade2 = as.integer(st$Tgrade == 2),
                  Grade3 = as.integer(st$Tgrade == 3),
                  TD = as.integer(st$TX_Delay == SHORT_CODE),
                  BD = as.integer(st$BX_Delay == SHORT_CODE))
  d[!is.na(d$county) & d$county >= 1 & d$county <= N & !is.na(d$time) &
    d$time > 0 & !is.na(d$delta), ]
}
## check that the rebuilt stratum is row-for-row the one analysed
check_aligned <- function(d, x) {
  ok <- nrow(d) == length(x$delta) && all(d$delta == x$delta) && all(d$county == x$county) &&
    isTRUE(all.equal(unname(as.matrix(d[, c("Age", "HR", "Grade2", "Grade3", "TD", "BD")])),
                     unname(x$Xraw), check.attributes = FALSE))
  if (!ok) stop("data do not match the saved outputs (different data file or SHORT_CODE?)")
}

## BRFSS counts, mapped to counties exactly as in fcr_data_analysis.R
brfss <- function(race) {
  bf <- suppressMessages(read_excel(if (race == "AA") "mammo_county_16A.xlsx"
                                    else "mammo_county_16W.xlsx", skip = 4, col_names = FALSE))
  bcode <- suppressWarnings(as.numeric(bf[[1]]))
  bcty  <- ifelse(!is.na(bcode) & bcode == 86, 13, (bcode + 1)/2)   # Miami-Dade: FIPS 086 -> county 13
  b_no  <- suppressWarnings(as.numeric(bf[[2]]))
  b_yes <- suppressWarnings(as.numeric(bf[[3]]))
  ok <- !is.na(bcty) & bcty >= 1 & bcty <= N
  m0 <- n0 <- rep(0, N)
  m0[round(bcty[ok])] <- b_yes[ok]; n0[round(bcty[ok])] <- b_no[ok] + b_yes[ok]
  list(m0 = m0, n0 = n0)
}
any_run <- function(race) { k <- grep(paste0("^", race, "_"), have, value = TRUE)
                            if (length(k)) runs[[k[1]]] else NULL }

################################################################################
# Section 5 text: stratified log-rank test for racial differences
################################################################################
lr <- survdiff(Surv(yr, delta) ~ race + strata(stage), data = surv_raw)
p_lr <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
sink(fp("logrank.txt"))
cat("Stratified (by stage) log-rank test, AA vs WA\n")
print(lr)
cat(sprintf("\np-value = %.3g\n\nBy stage:\n", p_lr))
for (s in levels(surv_raw$stage)) {
  l <- survdiff(Surv(yr, delta) ~ race, data = surv_raw[surv_raw$stage == s, ])
  cat(sprintf("  %-9s chisq = %8.2f  p = %.3g\n", s, l$chisq, pchisq(l$chisq, 1, lower.tail = FALSE)))
}
sink()

################################################################################
# Table 1: posterior variable inclusion proportions
################################################################################
covs <- c("Age", "HR", "Grade", "SMU", "TD", "BD")
tab1 <- sapply(STRATA, function(k) if (is.null(runs[[k]])) rep(NA, 6) else runs[[k]]$vip6_weighted[covs])
rownames(tab1) <- c("Age", "HR Status", "Tumor Grade", "SMU", "Treatment Delay", "Biopsy Delay")
colnames(tab1) <- c("(AA, L)", "(WA, L)", "(AA, R)", "(WA, R)", "(AA, D)", "(WA, D)")
write.csv(round(tab1, 3), fp("table1_vip.csv"))
tex_cell <- function(v, j) {                    # top two (after rounding) in bold
  if (is.na(v[j])) return("--")
  r <- round(v, 3); s <- f3(v[j])
  if (r[j] >= sort(r, decreasing = TRUE)[2]) sprintf("\\textbf{%s}", s) else s
}
tex1 <- sapply(seq_len(ncol(tab1)), function(c) sapply(seq_len(nrow(tab1)), function(r) tex_cell(tab1[, c], r)))
writeLines(c(paste0(paste(c("Strata: (Race, Stage)", colnames(tab1)), collapse = " & "), " \\\\ \\hline"),
             paste0(apply(cbind(rownames(tab1), tex1), 1, paste, collapse = " & "), " \\\\")),
           fp("table1_vip.tex"))
cat("Table 1 written\n"); missing_note(STRATA)

################################################################################
# Figure 3: county proportions of long BD / long TD and No SMU (1 - M_hat), by race
################################################################################
box_vals <- list()
for (r in c("AA", "WA")) {
  d <- surv_raw[surv_raw$race == r, ]
  box_vals[[paste("BD", r)]]  <- tapply(d$BX_Delay != SHORT_CODE, d$county, mean)
  box_vals[[paste("SMU", r)]] <- if (!is.null(any_run(r))) 1 - any_run(r)$M_hat else NA
  box_vals[[paste("TD", r)]]  <- tapply(d$TX_Delay != SHORT_CODE, d$county, mean)
}
box_vals <- box_vals[c("BD AA", "BD WA", "SMU AA", "SMU WA", "TD AA", "TD WA")]
pdf(fp("figure3_boxplots.pdf"), width = 7, height = 4.5)
par(mar = c(4.5, 4.2, 1, 1))
boxplot(box_vals, at = c(1, 2, 4, 5, 7, 8), names = rep(c("AA", "WA"), 3),
        ylab = "Proportion across counties", col = rep(c("grey70", "white"), 3))
axis(1, at = c(1.5, 4.5, 7.5), tick = FALSE, line = 1.5,
     labels = c("Long BD", expression("No SMU (1 - " * hat(M)[i] * ")"), "Long TD"))
invisible(dev.off())

################################################################################
# Figure 4: estimated survival of representative patients in Broward
################################################################################
pdf(fp("figure4_survival.pdf"), width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2, 1))
for (s in c("L", "R", "D")) {
  plot(NA, xlim = c(0, 10), ylim = c(0, 1), xlab = "Years", ylab = expression(hat(S)(t)),
       main = c(L = "Local", R = "Regional", D = "Distant")[s])
  for (r in c("AA", "WA")) { x <- runs[[paste0(r, "_", s)]]
    if (!is.null(x)) lines(c(0, x$times), c(1, x$S_rep[1, ]), lty = if (r == "AA") 2 else 1, lwd = 1.5) }
  legend("topright", c("AA", "WA"), lty = c(2, 1), bty = "n")
}
invisible(dev.off())

################################################################################
# Table 2 (DMS(10), AA Broward, Local and Regional) and Table S2 (DMS(5))
################################################################################
fmt_ci <- function(m, l, h) sprintf("%s (%s, %s)", f3(m), f3(l), f3(h))
lab     <- c(SMU = "SMU (M* = M_max vs M_hat)", BD = "Biopsy Delay (short vs long)",
             TD = "Treatment Delay (short vs long)")
lab_tex <- c(SMU = "SMU: $M_i^*=M^{\\max}$ vs $M_i=\\widehat{M}_i$", BD = "Biopsy Delay: short vs long",
             TD = "Treatment Delay: short vs long")
tab2_stage <- function(k) { x <- runs[[k]]; if (is.null(x)) return(rep(NA, 9))
  t2 <- x$tab2[order(match(x$tab2$contrast, names(lab)), x$tab2$grade), ]
  fmt_ci(t2$DMS10, t2$HPD10_lo, t2$HPD10_hi) }
tab2 <- data.frame(Covariate = rep(lab, each = 3), Grade = rep(1:3, 3),
                   Local = tab2_stage("AA_L"), Regional = tab2_stage("AA_R"))
write.csv(tab2, fp("table2_dms10.csv"), row.names = FALSE)
writeLines(c("Covariate & Tumor Grade & Local & Regional \\\\ \\hline",
             paste0(apply(cbind(rep(lab_tex, each = 3), tab2[, -1]), 1, paste, collapse = " & "), " \\\\")),
           fp("table2_dms10.tex"))
s2_stage <- function(k) { x <- runs[[k]]; if (is.null(x)) return(rep(NA, 6))
  t2 <- x$tab2[x$tab2$contrast %in% c("BD", "SMU"), ]
  t2 <- t2[order(match(t2$contrast, c("BD", "SMU")), t2$grade), ]
  f3(t2$DMS5) }
tabS2 <- data.frame(Covariate = rep(c("Biopsy Delay", "SMU"), each = 3), Grade = rep(1:3, 2),
                    Local = s2_stage("AA_L"), Regional = s2_stage("AA_R"), Distant = s2_stage("AA_D"))
write.csv(tabS2, fp("tableS2_dms5.csv"), row.names = FALSE)
writeLines(c("Covariate & Tumor Grade & Local & Regional & Distant \\\\ \\hline",
             paste0(apply(tabS2, 1, paste, collapse = " & "), " \\\\")), fp("tableS2_dms5.tex"))
cat("Tables 2 and S2 written\n"); missing_note(c("AA_L", "AA_R", "AA_D"))

################################################################################
# Per-county DMS(10), AA Regional stage: Figures 5, S9-S12
################################################################################
ct <- if (!is.null(runs$AA_R)) runs$AA_R$cty_tab else NULL
if (!is.null(ct)) {
  ct$name <- COUNTY[ct$county]
  write.csv(ct[, c("county", "name", "n_i", "M_hat", "SMU10", "SMU10_hpd_lo", "SMU10_hpd_hi",
                   "TD10", "TD10_hpd_lo", "TD10_hpd_hi", "BD10", "BD10_hpd_lo", "BD10_hpd_hi")],
            fp("county_dms10_AA_R.csv"), row.names = FALSE)

  county_plot <- function(file, y, ylab) {
    pdf(fp(file), width = 7.2, height = 5); par(mar = c(4.2, 4.6, 1, 1), las = 1)
    hl <- ct$county %in% HIGHLIGHT
    plot(ct$n_i[!hl], y[!hl], pch = 16, cex = 0.6, xlab = expression(n[i]), ylab = ylab,
         xlim = range(ct$n_i), ylim = range(y))
    abline(h = 0, lty = 3)
    for (nm in names(HIGHLIGHT)) { j <- ct$county == HIGHLIGHT[nm]
      if (any(j)) points(ct$n_i[j], y[j], pch = HL_PCH[nm], cex = 1.1, lwd = 1.2) }
    legend("topright", names(HIGHLIGHT), pch = HL_PCH[names(HIGHLIGHT)], bty = "n")
    invisible(dev.off())
  }
  county_plot("figure5_dms10_smu.pdf", ct$SMU10, expression(DMS(10) * ":  " * M[max] * " vs " * hat(M)[i]))
  county_plot("figS11_dms10_td.pdf", ct$TD10, "DMS(10): short vs long TD")

  ## Figure S12: posterior mean and 95% HPD interval against n_i,
  ## counties with n_i < 5 and with 50 < n_i < 100
  sel <- ct[ct$n_i < 5 | (ct$n_i > 50 & ct$n_i < 100), ]
  xs <- jitter(sel$n_i, amount = 0.4)
  w <- sel$TD10_hpd_hi > sel$TD10_hpd_lo
  pdf(fp("figS12_dms10_td_ci.pdf"), width = 7.2, height = 5); par(mar = c(4.2, 4.2, 1, 1), las = 1)
  plot(xs, sel$TD10, ylim = range(sel$TD10_hpd_lo, sel$TD10_hpd_hi), pch = 16,
       xlab = expression(n[i]), ylab = "DMS(10): short vs long TD")
  arrows(xs[w], sel$TD10_hpd_lo[w], xs[w], sel$TD10_hpd_hi[w], angle = 90, code = 3, length = 0.03)
  abline(h = 0, lty = 3)
  invisible(dev.off())

  ## Figures S9, S10: maps (package 'maps')
  if (requireNamespace("maps", quietly = TRUE)) {
    m <- maps::map("county", "florida", fill = TRUE, plot = FALSE)
    key <- function(s) gsub("[^a-z]", "", tolower(s))
    poly_cty <- match(key(sub(":.*", "", sub("^florida,", "", m$names))), key(COUNTY))
    if (anyNA(poly_cty)) warning("map polygons not matched: ", paste(m$names[is.na(poly_cty)], collapse = ", "))
    draw_map <- function(file, v, title) {
      val <- rep(NA_real_, N); val[ct$county] <- v
      brk <- pretty(range(val, na.rm = TRUE), 6)
      pal <- grey.colors(length(brk) - 1, start = 0.95, end = 0.25)
      col <- pal[cut(val[poly_cty], brk, include.lowest = TRUE)]; col[is.na(col)] <- "white"
      pdf(fp(file), width = 7, height = 6)
      maps::map("county", "florida", fill = TRUE, col = col, mar = c(1, 1, 2, 1))
      title(title)
      legend("bottomleft", legend = c(levels(cut(0, brk, include.lowest = TRUE)), "no AA Regional patients"),
             fill = c(pal, "white"), title = "DMS(10), years", bty = "n", cex = 0.8)
      invisible(dev.off())
    }
    draw_map("figS9_map_dms10_smu.pdf", ct$SMU10, "DMS(10): M_max vs M_hat, AA Regional")
    draw_map("figS10_map_dms10_td.pdf", ct$TD10, "DMS(10): short vs long TD, AA Regional")
  } else cat("  Figures S9-S10 skipped: install.packages('maps') to draw the maps\n")
  cat("Figures 5, S9-S12 written\n")
} else missing_note("AA_R")

################################################################################
# Figure S5: Kaplan-Meier by stage and race
################################################################################
pdf(fp("figS5_km.pdf"), width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2, 1))
for (s in levels(surv_raw$stage)) {
  km <- survfit(Surv(yr, delta) ~ race, data = surv_raw[surv_raw$stage == s, ])
  plot(km, lty = c(2, 1), xlab = "Years", ylab = expression(hat(S)[KM](t)), main = s)
  legend("topright", c("AA", "WA"), lty = c(2, 1), bty = "n")
}
invisible(dev.off())

################################################################################
# Figures S6-S8: residuals, spatial SBART vs PH-frailty (Regional stage)
# PH-frailty comparator: Cox model with the Step-1 M_hat as a covariate and
# independent Gaussian county frailties, fitted by penalised partial likelihood
# (survival::coxph).
################################################################################
fit_ph <- function(d, M_hat) {
  d$Mhat <- M_hat[d$county]
  coxph(Surv(time, delta) ~ Age + HR + Grade2 + Grade3 + TD + BD + Mhat +
          frailty(county, distribution = "gaussian"), data = d,
        control = coxph.control(iter.max = 1000))
}
resid <- list()
for (k in c("AA_R", "WA_R")) { x <- runs[[k]]; if (is.null(x)) next
  d <- stratum_data(substr(k, 1, 2), substr(k, 4, 4))
  check_aligned(d, x)
  ph <- fit_ph(d, x$M_hat)
  resid[[k]] <- list(d = d, cs_sb = x$cox_snell, dv_sb = x$deviance_res,
                     cs_ph = d$delta - residuals(ph, type = "martingale"),
                     dv_ph = residuals(ph, type = "deviance"))
}
if (length(resid)) {
  pdf(fp("figS6_coxsnell.pdf"), width = 9, height = 4.5); par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  for (k in names(resid)) { r <- resid[[k]]
    na_sb <- survfit(Surv(r$cs_sb, r$d$delta) ~ 1); na_ph <- survfit(Surv(r$cs_ph, r$d$delta) ~ 1)
    lim <- c(0, max(na_sb$time, na_ph$time, na_sb$cumhaz, na_ph$cumhaz))
    plot(na_sb$time, na_sb$cumhaz, type = "s", col = "darkgreen", xlim = lim, ylim = lim,
         xlab = "Cox-Snell residual", ylab = "Nelson-Aalen cumulative hazard",
         main = c(AA_R = "AA, Regional", WA_R = "WA, Regional")[k])
    lines(na_ph$time, na_ph$cumhaz, type = "s", col = "red"); abline(0, 1, lty = 2)
    legend("topleft", c("Spatial SBART", "PH-frailty"), col = c("darkgreen", "red"), lty = 1, bty = "n")
  }
  invisible(dev.off())

  share <- function(v, ev) c(all = mean(abs(v) > 2), uncensored = mean(abs(v[ev]) > 2))
  shares <- do.call(rbind, lapply(names(resid), function(k) { r <- resid[[k]]; ev <- r$d$delta == 1
    data.frame(stratum = k, method = c("Spatial SBART", "PH-frailty"),
               rbind(share(r$dv_sb, ev), share(r$dv_ph, ev)),
               above_2 = c(sum(r$dv_sb > 2), sum(r$dv_ph > 2)), below_m2 = c(sum(r$dv_sb < -2), sum(r$dv_ph < -2))) }))
  write.csv(shares, fp("residual_shares.csv"), row.names = FALSE)

  for (k in names(resid)) { r <- resid[[k]]
    yl <- range(-2.5, 2.5, r$dv_sb, r$dv_ph)            # one scale for both methods
    pdf(fp(sprintf("figS%d_deviance_%s.pdf", if (k == "AA_R") 7 else 8, k)), width = 15, height = 10.5)
    par(mfrow = c(3, 4), mar = c(4.2, 4.2, 2.5, 1))
    for (v in c("Age", "HR", "Grade2", "Grade3", "TD", "BD")) {
      for (m in c("sb", "ph")) {
        dv <- r[[paste0("dv_", m)]]
        ttl <- sprintf("%s: %.1f%% outside (-2, 2)", c(sb = "Spatial SBART", ph = "PH-frailty")[m], 100*mean(abs(dv) > 2))
        if (v == "Age") plot(r$d$Age, dv, pch = 16, cex = 0.4, ylim = yl, xlab = "Age (years)",
                             ylab = "Deviance residual", main = ttl)
        else boxplot(dv ~ factor(r$d[[v]], levels = 0:1), ylim = yl, xlab = v,
                     ylab = "Deviance residual", main = ttl)
        abline(h = c(-2, 2), col = "red", lty = 2)
      }
    }
    invisible(dev.off())
  }
  cat("Figures S6-S8 written\n")
}
missing_note(c("AA_R", "WA_R"))

################################################################################
# Figures S13-S15: Step 1 for M (AA BRFSS); CAR model from fcr_data_analysis.R,
# spatially independent model (iid normal prior on logit M) fitted here
################################################################################
step1_iid <- function(m0, n0, n_iter = 12000, burn = 6000, seed = 20260816) {
  set.seed(seed); has0 <- n0 > 0
  u <- qlogis((m0 + 0.5)/(n0 + 1)); s0 <- 1; keep <- matrix(NA_real_, n_iter, length(u))
  for (it in seq_len(n_iter)) {
    for (i in seq_along(u)) {
      ui <- rnorm(1, u[i], 0.3)
      ll <- if (has0[i]) dbinom(m0[i], n0[i], plogis(ui), log = TRUE) -
                         dbinom(m0[i], n0[i], plogis(u[i]), log = TRUE) else 0
      if (log(runif(1)) < ll - 0.5*(ui^2 - u[i]^2)/s0) u[i] <- ui
    }
    s0 <- 1/rgamma(1, 1 + length(u)/2, 1 + sum(u^2)/2)
    keep[it, ] <- u
  }
  keep <- keep[(burn + 1):n_iter, , drop = FALSE]
  list(M_hat = colMeans(plogis(keep)), var_logitM = apply(keep, 2, var))
}
car <- if (!is.null(runs$AA_R)) runs$AA_R else any_run("AA")
if (!is.null(car)) {
  b <- brfss("AA"); ok <- b$n0 > 0
  iid <- step1_iid(b$m0, b$n0)
  pdf(fp("figS13_var_logitM.pdf"), width = 6, height = 5); par(mar = c(4.2, 4.6, 1, 1))
  plot(log(b$n0[ok]), car$var_logitM[ok], pch = 1, ylim = range(car$var_logitM[ok], iid$var_logitM[ok]),
       xlab = expression(log(n[0 * i])), ylab = expression(widehat(Var) * "{logit(" * M[i] * ") | " * D[0] * "}"))
  points(log(b$n0[ok]), iid$var_logitM[ok], pch = 2, col = "red")
  legend("topright", c("Spatial CAR", "Independent"), pch = 1:2, col = c("black", "red"), bty = "n")
  invisible(dev.off())
  pdf(fp("figS14_Mhat_comparison.pdf"), width = 5.5, height = 5.5); par(mar = c(4.2, 4.2, 1, 1))
  plot(car$M_hat, iid$M_hat, xlab = expression(hat(M)[i] * ", spatial CAR model"),
       ylab = expression(hat(M)[i] * ", independent model"))
  abline(0, 1, lty = 2)
  invisible(dev.off())
  pdf(fp("figS15_Mhat_ci.pdf"), width = 6, height = 5); par(mar = c(4.2, 4.2, 1, 1))
  q <- car$M_post_q; w <- ok & q[3, ] > q[1, ]
  plot(log(b$n0[ok]), car$M_hat[ok], pch = 16, ylim = range(q[, ok]),
       xlab = expression(log(n[0 * i])), ylab = expression(hat(M)[i] ~ "and 95% interval"))
  arrows(log(b$n0[w]), q[1, w], log(b$n0[w]), q[3, w], angle = 90, code = 3, length = 0.02, col = "grey40")
  invisible(dev.off())
  cat("Figures S13-S15 written", if (any(!ok)) paste0("(counties with n_0i = 0, not plotted on the log scale: ",
      paste(COUNTY[!ok], collapse = ", "), ")"), "\n")
} else missing_note("AA_R")

################################################################################
# Figure S16: frailty W_i, posterior mean and 95% interval, AA Regional
# (three counties: one with n_i < 5 and the smallest and largest with 50 < n_i < 100)
################################################################################
if (!is.null(runs$AA_R)) {
  W <- runs$AA_R$W_tab
  small <- which(W$n_i > 0 & W$n_i < 5); mid <- which(W$n_i > 50 & W$n_i < 100)
  idx <- c(small[which.max(W$n_i[small])], mid[which.min(W$n_i[mid])], mid[which.max(W$n_i[mid])])
  pick <- W[unique(idx[!is.na(idx)]), ]
  w <- pick$W_q975 > pick$W_q025
  pdf(fp("figS16_W_ci.pdf"), width = 6, height = 5); par(mar = c(4.2, 4.2, 1, 1))
  plot(pick$n_i, pick$W_mean, pch = 16, xlim = c(0, max(pick$n_i) * 1.08),
       ylim = range(pick$W_q025, pick$W_q975) * c(1, 1.1),
       xlab = expression(n[i]), ylab = expression(hat(W)[i] ~ "and 95% interval"))
  arrows(pick$n_i[w], pick$W_q025[w], pick$n_i[w], pick$W_q975[w], angle = 90, code = 3, length = 0.05)
  text(pick$n_i, pick$W_q975, COUNTY[pick$county], pos = 3, cex = 0.8)
  invisible(dev.off())
} else missing_note("AA_R")

################################################################################
# Figure S17: trace plots, AA Distant, Step 2 (chains in out/AA_D and out/AA_D_chain2)
################################################################################
chains <- Filter(Negate(is.null), list(runs$AA_D, load_run("AA_D", "AA_D_chain2")))
if (length(chains) == 2) {
  s1 <- chains[[1]]$settings; s2 <- chains[[2]]$settings
  if (s1["num_iter"] != s2["num_iter"] || s1["burn_in"] != s2["burn_in"])
    stop("the two AA_D chains differ in num_iter / burn_in")
  if (!is.na(s1["seed"]) && !is.na(s2["seed"]) && s1["seed"] == s2["seed"])
    warning("both AA_D chains used the same seed, so their traces are identical")
  if (chains[[2]]$settings[["use_step3"]] != 1) stop("AA_D_chain2 was run without Step 3")
  if (!isTRUE(all.equal(chains[[1]]$M_hat, chains[[2]]$M_hat)))
    stop("the two AA_D chains do not share the Step-1 output (different n1/b1 or data?)")
}
if (length(chains)) {
  x1 <- chains[[1]]; post <- (x1$settings["burn_in"] + 1):x1$settings["num_iter"]
  it <- post - x1$settings["burn_in"]
  pdf(fp("figS17_traces.pdf"), width = 12, height = 6.5); par(mfrow = c(2, 3), mar = c(4.2, 4.6, 2, 1))
  panel <- function(ys, ylab, main, tag) {
    plot(it, ys[[1]], type = "l", col = "blue", ylim = range(unlist(ys), na.rm = TRUE),
         xlab = "Iteration (post burn-in)", ylab = ylab, main = main)
    if (length(ys) > 1) lines(it, ys[[2]], col = "red")
    mtext(sprintf("(%s)", tag), side = 3, adj = 0, line = 0.3, font = 2)
  }
  panel(lapply(chains, function(x) log(x$lam_trace[post])), expression(log(lambda[0])), "", "a")
  panel(lapply(chains, function(x) log(x$s1_trace[post])), expression(log(sigma[1]^2)), "", "b")
  panel(lapply(chains, function(x) x$rho1_trace[post]), expression(rho[1]), "", "c")
  for (j in seq_along(TRACE_COUNTIES)) { c3 <- TRACE_COUNTIES[j]
    panel(lapply(chains, function(x) x$R_trace[post, c3]), bquote(log(W[.(c3)])), COUNTY[c3], letters[3 + j]) }
  invisible(dev.off())
  if (length(chains) < 2 && !is.null(runs$AA_D))
    cat("  Figure S17 drawn with one chain; run AA D again with another seed into",
        file.path(outdir, "AA_D_chain2"), "for the second chain\n")
} else missing_note("AA_D")

cat("\nAll outputs in", figdir, "\n")
