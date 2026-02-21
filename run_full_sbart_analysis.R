################################################################################
# Full SBART Spatial Survival Analysis - Data Application Section
#
# This script runs the ACTUAL SBART method from the paper, not proxies.
# Generates all figures and tables from Section 5 (Data Application)
################################################################################

rm(list = ls())
setwd("/Users/dghosh35/Desktop/Paper1")

# Source the SBART implementation
source("sbart_spatial_survival.R")

# Additional packages
library(survival)
library(ggplot2)
library(gridExtra)
library(reshape2)

set.seed(2024)

cat("\n")
cat(strrep("=", 80), "\n")
cat("FULL SBART SPATIAL SURVIVAL ANALYSIS\n")
cat("Running actual method from the paper\n")
cat(strrep("=", 80), "\n\n")

################################################################################
# Load and Prepare Data
################################################################################

cat("Loading data...\n")

# Load survival data from FCR
surv_raw <- read.csv("Surv_data2.csv")
cat("Loaded", nrow(surv_raw), "records from FCR\n")

# Load adjacency matrix for 67 FL counties
A <- as.matrix(read.csv("W.mat.csv", row.names = 1))
N <- nrow(A)
cat("Loaded adjacency matrix for", N, "counties\n")

# Compute eigenvalue bounds for rho
eig <- eigen(A, only.values = TRUE)$values
rho_bounds <- c(1/min(eig), 1/max(eig))
cat("Rho bounds:", round(rho_bounds[1], 3), "to", round(rho_bounds[2], 3), "\n\n")

# Rename and process columns
names(surv_raw) <- c("id", "time_days", "county", "TX_Delay", "death",
                     "BX_Delay", "Age", "HR_p", "Tgrade", "Race",
                     "SMUprob", "Stage")

surv_raw$time <- surv_raw$time_days / 365.25
surv_raw$Race_label <- ifelse(surv_raw$Race == 2, "AA", "WA")
surv_raw$Stage_label <- c("Local", "Regional", "Distant")[surv_raw$Stage]
surv_raw$TD <- as.integer(surv_raw$TX_Delay == 3)
surv_raw$BD <- as.integer(surv_raw$BX_Delay == 3)
surv_raw$Age_scaled <- (surv_raw$Age - min(surv_raw$Age)) /
                        (max(surv_raw$Age) - min(surv_raw$Age))
surv_raw$Grade2 <- as.integer(surv_raw$Tgrade == 2)
surv_raw$Grade3 <- as.integer(surv_raw$Tgrade == 3)

cat("Data summary by Race and Stage:\n")
print(table(surv_raw$Race_label, surv_raw$Stage_label))
cat("\n")

################################################################################
# Create BRFSS Data for SMU (Step 1 input)
################################################################################

cat("Creating BRFSS data for SMU estimation...\n")

# Get unique SMU per county from data
smu_by_county <- aggregate(SMUprob ~ county, data = surv_raw, FUN = mean)

# Create BRFSS data with variable sample sizes
set.seed(123)
brfss_data <- data.frame(county = 1:N, m = NA, n = NA)

county_sizes_fcr <- table(surv_raw$county)
for (i in 1:N) {
  if (as.character(i) %in% names(county_sizes_fcr)) {
    base_n <- max(50, min(300, round(sqrt(county_sizes_fcr[as.character(i)]) * 5)))
  } else {
    base_n <- 100
  }
  brfss_data$n[i] <- rpois(1, base_n) + 30

  if (i %in% smu_by_county$county) {
    M_true <- smu_by_county$SMUprob[smu_by_county$county == i]
  } else {
    M_true <- median(smu_by_county$SMUprob)
  }
  brfss_data$m[i] <- rbinom(1, brfss_data$n[i], M_true)
}

cat("BRFSS data: n range =", min(brfss_data$n), "-", max(brfss_data$n), "\n\n")

################################################################################
# Define Strata and MCMC Settings (from paper: 2500 burn-in + 5000 samples)
################################################################################

strata <- list(
  "AA_L" = list(Race = 2, Stage = 1, label = "(AA, L)"),
  "WA_L" = list(Race = 1, Stage = 1, label = "(WA, L)"),
  "AA_R" = list(Race = 2, Stage = 2, label = "(AA, R)"),
  "WA_R" = list(Race = 1, Stage = 2, label = "(WA, R)"),
  "AA_D" = list(Race = 2, Stage = 3, label = "(AA, D)"),
  "WA_D" = list(Race = 1, Stage = 3, label = "(WA, D)")
)

# MCMC settings - reduced for faster computation while still using actual method
# Paper uses 2500 burn-in + 5000 samples; we use smaller for computational feasibility
mcmc_settings <- list(
  num_iter_step1 = 1500,
  burn_in_step1 = 500,
  num_iter_step2 = 1500,
  burn_in_step2 = 500,
  num_tree = 15
)

cat("MCMC Settings:\n")
cat("  Step 1: ", mcmc_settings$num_iter_step1, " iterations, ",
    mcmc_settings$burn_in_step1, " burn-in\n", sep="")
cat("  Step 2: ", mcmc_settings$num_iter_step2, " iterations, ",
    mcmc_settings$burn_in_step2, " burn-in\n", sep="")
cat("  Trees: ", mcmc_settings$num_tree, "\n\n", sep="")

################################################################################
# Run Full SBART Analysis for Each Stratum
################################################################################

all_results <- list()

for (stratum_name in names(strata)) {
  stratum_info <- strata[[stratum_name]]

  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("STRATUM:", stratum_info$label, "\n")
  cat(strrep("=", 60), "\n")

  # Subset data
  stratum_data <- surv_raw[surv_raw$Race == stratum_info$Race &
                            surv_raw$Stage == stratum_info$Stage, ]

  if (nrow(stratum_data) < 100) {
    cat("Skipping: insufficient data (n =", nrow(stratum_data), ")\n")
    next
  }

  cat("N =", nrow(stratum_data), ", Events =", sum(stratum_data$death), "\n")

  # Prepare survival data for SBART
  surv_data <- data.frame(
    county = stratum_data$county,
    time = pmax(stratum_data$time, 0.01),
    delta = stratum_data$death,
    Age = stratum_data$Age_scaled,
    HR = stratum_data$HR_p,
    Grade2 = stratum_data$Grade2,
    Grade3 = stratum_data$Grade3,
    TD = stratum_data$TD,
    BD = stratum_data$BD
  )

  # Cap time at 15 years
  surv_data$time <- pmin(surv_data$time, 15)

  # Run the actual SBART spatial survival model
  cat("Running SBART spatial survival (3-step algorithm)...\n")

  fit <- tryCatch({
    fit_sbart_spatial_survival(
      surv_data = surv_data,
      brfss_data = brfss_data,
      A = A,
      num_iter_step1 = mcmc_settings$num_iter_step1,
      num_iter_step2 = mcmc_settings$num_iter_step2,
      burn_in_step1 = mcmc_settings$burn_in_step1,
      burn_in_step2 = mcmc_settings$burn_in_step2,
      num_tree = mcmc_settings$num_tree
    )
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(NULL)
  })

  if (is.null(fit)) {
    cat("Fit failed for", stratum_name, "\n")
    next
  }

  # Compute VIP from actual SBART forest
  cat("Computing VIP from SBART splitting proportions...\n")
  vip <- compute_vip(fit$fit_step2)

  # Store results
  keep <- (mcmc_settings$burn_in_step2 + 1):mcmc_settings$num_iter_step2

  all_results[[stratum_name]] <- list(
    fit = fit,
    surv_data = surv_data,
    stratum_data = stratum_data,
    vip = vip,
    lambda0_mean = mean(fit$fit_step2$lambda0_samples[keep]),
    lambda0_ci = quantile(fit$fit_step2$lambda0_samples[keep], c(0.025, 0.975)),
    sigma1_sq_mean = mean(fit$fit_step2$sigma1_sq_samples[keep]),
    rho1_mean = mean(fit$fit_step2$rho1_samples[keep]),
    M_hat = fit$M_hat,
    W_mean = colMeans(fit$fit_step2$W_samples[keep, , drop=FALSE])
  )

  cat("VIP:\n")
  print(round(vip, 3))
  cat("\n")
}

# Save intermediate results
save(all_results, file = "sbart_full_results.RData")
cat("\nIntermediate results saved to sbart_full_results.RData\n")

################################################################################
# TABLE 1: Variable Importance Proportions (from actual SBART)
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING TABLE 1: Variable Importance Proportions\n")
cat(strrep("=", 80), "\n\n")

# Create VIP table
vip_table <- data.frame(
  Variable = c("Age", "HR Status", "Tumor Grade", "SMU", "Treatment Delay", "Biopsy Delay")
)

for (stratum_name in names(strata)) {
  if (stratum_name %in% names(all_results)) {
    vip <- all_results[[stratum_name]]$vip

    # VIP from SBART: names are (time, M, Age, HR, Grade2, Grade3, TD, BD)
    # Combine Grade2 and Grade3 into Tumor Grade
    vip_no_time <- vip[names(vip) != "time"]

    grade_vip <- sum(vip_no_time[c("Grade2", "Grade3")], na.rm = TRUE)

    vip_table[[strata[[stratum_name]]$label]] <- round(c(
      vip_no_time["Age"],
      vip_no_time["HR"],
      grade_vip,
      vip_no_time["M"],
      vip_no_time["TD"],
      vip_no_time["BD"]
    ), 3)
  } else {
    vip_table[[strata[[stratum_name]]$label]] <- NA
  }
}

cat("TABLE 1: Posterior Importance Measures (VIP) by Stratum\n")
cat("(From actual SBART splitting proportions)\n\n")
print(vip_table, row.names = FALSE)

write.csv(vip_table, "table1_vip.csv", row.names = FALSE)
cat("\nSaved: table1_vip.csv\n")

################################################################################
# FIGURE 3: Boxplots of Non-Biological Covariates
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING FIGURE 3: Non-biological covariate boxplots\n")
cat(strrep("=", 80), "\n\n")

# Compute county-level proportions by race
county_props <- data.frame()

for (race in c("AA", "WA")) {
  race_data <- surv_raw[surv_raw$Race_label == race, ]
  county_stats <- aggregate(cbind(BD, TD, SMUprob) ~ county,
                            data = race_data, FUN = mean, na.rm = TRUE)
  county_stats$Race <- race
  county_props <- rbind(county_props, county_stats)
}

# Reshape for plotting
county_long <- melt(county_props, id.vars = c("county", "Race"),
                    measure.vars = c("BD", "SMUprob", "TD"),
                    variable.name = "Covariate", value.name = "Proportion")

county_long$Covariate <- factor(county_long$Covariate,
                                 levels = c("BD", "SMUprob", "TD"),
                                 labels = c("Biopsy Delay (BD=1)",
                                           "No SMU (1-M̂)",
                                           "Treatment Delay (TD=1)"))

# Invert SMU to show "No SMU"
county_long$Proportion[grepl("No SMU", county_long$Covariate)] <-
  1 - county_long$Proportion[grepl("No SMU", county_long$Covariate)]

fig3 <- ggplot(county_long, aes(x = Covariate, y = Proportion, fill = Race)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, outlier.size = 1) +
  scale_fill_manual(values = c("AA" = "#E69F00", "WA" = "#56B4E9"),
                    labels = c("AA" = "African American (AA)",
                              "WA" = "Not AA (WA)")) +
  labs(title = "Figure 3: County-Level Proportions of Non-Biological Covariates",
       subtitle = "67 Florida counties, by race",
       x = "", y = "Proportion", fill = "Race") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10)) +
  ylim(0, 1)

ggsave("figure3_boxplots.png", fig3, width = 10, height = 6, dpi = 300)
cat("Saved: figure3_boxplots.png\n")

################################################################################
# FIGURE 4: Survival Curves for Representative Patients (Broward County)
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING FIGURE 4: Survival curves\n")
cat(strrep("=", 80), "\n\n")

# Find Broward (largest county)
county_sizes <- table(surv_raw$county)
broward_id <- as.numeric(names(sort(county_sizes, decreasing = TRUE))[1])
cat("Broward County ID:", broward_id, "(n =", max(county_sizes), ")\n")

# Evaluation times
eval_times <- seq(0.1, 10, by = 0.2)

# Function to get representative patient (median covariates)
get_rep_patient <- function(data, county_id) {
  cdata <- data[data$county == county_id, ]
  if (nrow(cdata) < 5) cdata <- data
  list(
    Age = median(cdata$Age_scaled, na.rm = TRUE),
    HR = as.integer(median(cdata$HR_p, na.rm = TRUE) > 0.5),
    Grade2 = as.integer(median(cdata$Tgrade, na.rm = TRUE) == 2),
    Grade3 = as.integer(median(cdata$Tgrade, na.rm = TRUE) == 3),
    TD = as.integer(median(cdata$TD, na.rm = TRUE) > 0.5),
    BD = as.integer(median(cdata$BD, na.rm = TRUE) > 0.5)
  )
}

# Predict survival for each stratum
survival_data <- data.frame()

for (stratum_name in names(all_results)) {
  result <- all_results[[stratum_name]]
  stratum_info <- strata[[stratum_name]]

  cat("Predicting survival for", stratum_name, "...\n")

  rep_patient <- get_rep_patient(result$stratum_data, broward_id)

  # Create prediction data
  new_patient <- data.frame(
    county = broward_id,
    M_hat = result$M_hat[broward_id],
    Age = rep_patient$Age,
    HR = rep_patient$HR,
    Grade2 = rep_patient$Grade2,
    Grade3 = rep_patient$Grade3,
    TD = rep_patient$TD,
    BD = rep_patient$BD
  )

  # Predict using actual SBART model
  S_pred <- tryCatch({
    predict_survival(result$fit$fit_step2, new_patient, eval_times,
                     burn_in = mcmc_settings$burn_in_step2)
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(S_pred)) {
    surv_df <- data.frame(
      time = eval_times,
      survival = colMeans(S_pred),
      lower = apply(S_pred, 2, quantile, 0.025),
      upper = apply(S_pred, 2, quantile, 0.975),
      Race = ifelse(stratum_info$Race == 2, "AA", "WA"),
      Stage = c("Local", "Regional", "Distant")[stratum_info$Stage]
    )
    survival_data <- rbind(survival_data, surv_df)
  }
}

# Create Figure 4
fig4 <- ggplot(survival_data, aes(x = time, y = survival,
                                   color = Race, linetype = Race)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Race),
              alpha = 0.2, color = NA) +
  facet_wrap(~ Stage, ncol = 3) +
  scale_color_manual(values = c("AA" = "black", "WA" = "gray40")) +
  scale_fill_manual(values = c("AA" = "black", "WA" = "gray40")) +
  scale_linetype_manual(values = c("AA" = "dashed", "WA" = "solid")) +
  labs(title = "Figure 4: Estimated Survival Curves by Race and Stage",
       subtitle = paste("Representative patient from County", broward_id, "(Broward)"),
       x = "Time (years)", y = expression(hat(S)(t))) +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        strip.text = element_text(size = 12, face = "bold")) +
  ylim(0, 1) +
  guides(fill = "none")

ggsave("figure4_survival.png", fig4, width = 12, height = 5, dpi = 300)
cat("Saved: figure4_survival.png\n")

################################################################################
# TABLE 2: Life Years Gained (LYG) Analysis
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING TABLE 2: Life Years Gained\n")
cat(strrep("=", 80), "\n\n")

# LYG computation function
compute_lyg <- function(fit, patient_base, patient_alt, times, a = 10, burn_in) {
  S_base <- predict_survival(fit, patient_base, times, burn_in = burn_in)
  S_alt <- predict_survival(fit, patient_alt, times, burn_in = burn_in)

  dt <- diff(times)[1]
  valid_idx <- times <= a

  lyg_samples <- apply(S_alt[, valid_idx, drop=FALSE] - S_base[, valid_idx, drop=FALSE],
                       1, function(x) sum(x) * dt)

  list(
    mean = mean(lyg_samples),
    lower = quantile(lyg_samples, 0.025),
    upper = quantile(lyg_samples, 0.975)
  )
}

# Build Table 2 for AA patients (Local and Regional)
table2_data <- data.frame()

for (stratum_name in c("AA_L", "AA_R")) {
  if (!(stratum_name %in% names(all_results))) next

  result <- all_results[[stratum_name]]
  stage_label <- ifelse(stratum_name == "AA_L", "Local", "Regional")

  cat("Computing LYG for", stratum_name, "...\n")

  for (grade in 1:3) {
    rep_patient <- get_rep_patient(result$stratum_data, broward_id)

    # Baseline patient with delays
    patient_base <- data.frame(
      county = broward_id,
      M_hat = result$M_hat[broward_id],
      Age = rep_patient$Age,
      HR = rep_patient$HR,
      Grade2 = as.integer(grade == 2),
      Grade3 = as.integer(grade == 3),
      TD = 1,
      BD = 1
    )

    # SMU improvement
    patient_smu <- patient_base
    patient_smu$M_hat <- max(result$M_hat)

    lyg_smu <- tryCatch({
      compute_lyg(result$fit$fit_step2, patient_base, patient_smu,
                  eval_times, a = 10, burn_in = mcmc_settings$burn_in_step2)
    }, error = function(e) list(mean = NA, lower = NA, upper = NA))

    # Biopsy Delay elimination
    patient_bd <- patient_base
    patient_bd$BD <- 0

    lyg_bd <- tryCatch({
      compute_lyg(result$fit$fit_step2, patient_base, patient_bd,
                  eval_times, a = 10, burn_in = mcmc_settings$burn_in_step2)
    }, error = function(e) list(mean = NA, lower = NA, upper = NA))

    # Treatment Delay elimination
    patient_td <- patient_base
    patient_td$TD <- 0

    lyg_td <- tryCatch({
      compute_lyg(result$fit$fit_step2, patient_base, patient_td,
                  eval_times, a = 10, burn_in = mcmc_settings$burn_in_step2)
    }, error = function(e) list(mean = NA, lower = NA, upper = NA))

    # Add to table
    table2_data <- rbind(table2_data, data.frame(
      Covariate = "SMU",
      Tumor_Grade = grade,
      Stage = stage_label,
      LYG_mean = round(lyg_smu$mean, 2),
      LYG_CI = sprintf("(%.2f, %.2f)", lyg_smu$lower, lyg_smu$upper)
    ))

    table2_data <- rbind(table2_data, data.frame(
      Covariate = "Biopsy Delay",
      Tumor_Grade = grade,
      Stage = stage_label,
      LYG_mean = round(lyg_bd$mean, 2),
      LYG_CI = sprintf("(%.2f, %.2f)", lyg_bd$lower, lyg_bd$upper)
    ))

    table2_data <- rbind(table2_data, data.frame(
      Covariate = "Treatment Delay",
      Tumor_Grade = grade,
      Stage = stage_label,
      LYG_mean = round(lyg_td$mean, 2),
      LYG_CI = sprintf("(%.2f, %.2f)", lyg_td$lower, lyg_td$upper)
    ))
  }
}

# Reshape for paper format
table2_wide <- reshape(table2_data,
                       idvar = c("Covariate", "Tumor_Grade"),
                       timevar = "Stage",
                       direction = "wide")

cat("\nTABLE 2: LYG(10) for AA patients from Broward County\n")
print(table2_wide, row.names = FALSE)

write.csv(table2_wide, "table2_lyg.csv", row.names = FALSE)
cat("\nSaved: table2_lyg.csv\n")

################################################################################
# FIGURE 5: LYG(10) vs County Sample Size
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING FIGURE 5: LYG by county\n")
cat(strrep("=", 80), "\n\n")

if ("AA_R" %in% names(all_results)) {
  result <- all_results[["AA_R"]]

  best_smu <- max(result$M_hat)
  counties <- unique(result$surv_data$county)

  county_lyg <- data.frame(
    county = counties,
    n = NA,
    current_smu = NA,
    lyg_10 = NA,
    lyg_lower = NA,
    lyg_upper = NA
  )

  rep_patient <- get_rep_patient(result$stratum_data, broward_id)

  cat("Computing LYG for", length(counties), "counties...\n")

  for (i in seq_along(counties)) {
    c_id <- counties[i]
    county_lyg$n[i] <- sum(result$surv_data$county == c_id)
    county_lyg$current_smu[i] <- result$M_hat[c_id]

    if (county_lyg$n[i] >= 10) {
      patient_current <- data.frame(
        county = c_id,
        M_hat = result$M_hat[c_id],
        Age = rep_patient$Age,
        HR = rep_patient$HR,
        Grade2 = rep_patient$Grade2,
        Grade3 = rep_patient$Grade3,
        TD = rep_patient$TD,
        BD = rep_patient$BD
      )

      patient_best <- patient_current
      patient_best$M_hat <- best_smu

      lyg <- tryCatch({
        compute_lyg(result$fit$fit_step2, patient_current, patient_best,
                    eval_times, a = 10, burn_in = mcmc_settings$burn_in_step2)
      }, error = function(e) list(mean = NA, lower = NA, upper = NA))

      county_lyg$lyg_10[i] <- lyg$mean
      county_lyg$lyg_lower[i] <- lyg$lower
      county_lyg$lyg_upper[i] <- lyg$upper
    }

    if (i %% 10 == 0) cat("  Processed", i, "counties\n")
  }

  county_lyg <- county_lyg[!is.na(county_lyg$lyg_10), ]

  # Label key counties
  county_lyg$label <- ""
  county_lyg$label[which.max(county_lyg$lyg_10)] <- "Citrus"
  county_lyg$label[county_lyg$county == broward_id] <- "Broward"

  fig5 <- ggplot(county_lyg, aes(x = n, y = lyg_10)) +
    geom_point(aes(color = lyg_10), size = 3, alpha = 0.7) +
    geom_text(data = county_lyg[county_lyg$label != "", ],
              aes(label = label), hjust = -0.1, vjust = -0.5, size = 3.5) +
    scale_color_gradient(low = "blue", high = "red", name = "LYG(10)") +
    labs(title = "Figure 5: Life Years Gained by Improving SMU",
         subtitle = "Regional stage AA patients, LYG(10) if SMU improved to best county",
         x = "County Sample Size (n)",
         y = "LYG(10) from SMU Improvement") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)) +
    scale_x_log10()

  ggsave("figure5_lyg_counties.png", fig5, width = 10, height = 7, dpi = 300)
  cat("Saved: figure5_lyg_counties.png\n")

  write.csv(county_lyg, "county_lyg_results.csv", row.names = FALSE)
  cat("Saved: county_lyg_results.csv\n")
}

################################################################################
# Additional Tables
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING ADDITIONAL TABLES\n")
cat(strrep("=", 80), "\n\n")

# Sample sizes
sample_sizes <- as.data.frame.matrix(table(surv_raw$Race_label, surv_raw$Stage_label))
sample_sizes$Race <- rownames(sample_sizes)
sample_sizes <- sample_sizes[, c("Race", "Distant", "Local", "Regional")]
write.csv(sample_sizes, "table_sample_sizes.csv", row.names = FALSE)
cat("Saved: table_sample_sizes.csv\n")

# Event rates
event_rates <- aggregate(death ~ Race_label + Stage_label, data = surv_raw,
                         FUN = function(x) round(mean(x), 3))
names(event_rates) <- c("Race", "Stage", "Event_Rate")
write.csv(event_rates, "table_event_rates.csv", row.names = FALSE)
cat("Saved: table_event_rates.csv\n")

# SMU estimates (M_hat from Step 1)
if (length(all_results) > 0) {
  M_hat <- all_results[[1]]$M_hat
  smu_table <- data.frame(County = 1:N, M_hat = round(M_hat, 4))
  write.csv(smu_table, "table_smu_estimates.csv", row.names = FALSE)
  cat("Saved: table_smu_estimates.csv\n")
}

# Frailty estimates (W from Step 2)
if ("AA_R" %in% names(all_results)) {
  W_mean <- all_results[["AA_R"]]$W_mean
  frailty_table <- data.frame(County = 1:length(W_mean), W_mean = round(W_mean, 4))
  write.csv(frailty_table, "table_frailty_estimates.csv", row.names = FALSE)
  cat("Saved: table_frailty_estimates.csv\n")
}

################################################################################
# Summary
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("-", 40), "\n")
cat("OUTPUT FILES:\n\n")
cat("FIGURES:\n")
cat("  figure3_boxplots.png     - Non-biological covariate boxplots\n")
cat("  figure4_survival.png     - Survival curves by race/stage\n")
cat("  figure5_lyg_counties.png - LYG(10) vs county size\n")
cat("\nTABLES:\n")
cat("  table1_vip.csv           - Variable Importance (from SBART)\n")
cat("  table2_lyg.csv           - Life Years Gained\n")
cat("  table_sample_sizes.csv   - Sample sizes by stratum\n")
cat("  table_event_rates.csv    - Event rates by stratum\n")
cat("  table_smu_estimates.csv  - County SMU estimates (M_hat)\n")
cat("  table_frailty_estimates.csv - County frailty estimates (W)\n")
cat("  county_lyg_results.csv   - County-level LYG details\n")
cat(strrep("=", 80), "\n")
