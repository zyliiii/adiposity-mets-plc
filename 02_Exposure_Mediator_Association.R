# ==============================================================================
# Script Name: 02_Exposure_Mediator_Association.R
# Purpose:     Linear regression for Exposure-Mediator (E-M) association
#              Main analysis and 3 sensitivity analyses included
# ==============================================================================

# Load required packages
library(rms)
library(dplyr)
library(purrr)

# Define variables
exp_vars <- c("zBMI", "zWaist", "zWHR", "zWHtR", "zABSI", "zHip", "zweightdiff")
met_vars <- paste0("zlog2met", 1:186)
covariates <- "startage + edu + inc + DIET + eversmoker + everdrinker + MET_total_wk"

# 1. Core Function: Run Linear and RCS Models for E-M Associations####
run_EM_analysis <- function(data, exposures, outcomes, covs) {
  # Update global datadist for the current dataset
  dd <<- datadist(data)
  options(datadist = "dd")
  # Iterate over each exposure
  results_list <- purrr::map(exposures, function(expo) {
    # Iterate over each metabolite outcome
    res_df <- purrr::map_dfr(outcomes, function(outc) {
      # (1) Linear Model Stats
      form_lin <- as.formula(paste(outc, "~", expo, "+", covs))
      fit_lin <- ols(form_lin, data = data)
      
      beta <- coef(fit_lin)[expo]
      se <- fit_lin$var[expo, expo]^0.5 
      t_val <- beta / se
      p_beta <- 2 * (1 - pt(abs(t_val), fit_lin$df.residual))
      
      beta_LCI <- beta - 1.96 * se
      beta_UCI <- beta + 1.96 * se
      beta_ci <- sprintf("%.2f (%.2f, %.2f)", beta, beta_LCI, beta_UCI)
      
      # (2) RCS Stats (Find best knots 3-5 based on min AIC)
      min_AIC <- Inf
      best_knots <- NA
      p_overall <- NA
      p_nonlin <- NA
      
      for (k in 3:5) {
        form_rcs <- as.formula(paste(outc, "~ rcs(", expo, ",", k, ") +", covs))
        fit_rcs <- ols(form_rcs, data = data)
        aic_val <- AIC(fit_rcs)
        
        if (aic_val < min_AIC) {
          min_AIC <- aic_val
          best_knots <- k
          anova_res <- anova(fit_rcs)
          p_overall <- anova_res[expo, "P"]  
          p_nonlin <- anova_res[" Nonlinear", "P"]
        }
      }
      # Return single row dataframe for this metabolite
      data.frame(
        metabolite = outc,
        EM_beta = beta,
        EM_beta_LCI = beta_LCI,
        EM_beta_UCI = beta_UCI,
        EM_beta_ci = beta_ci,
        EM_beta_P = p_beta,
        EM_best_knots = best_knots,
        EM_P_overall = p_overall,
        EM_P_nonlinear = p_nonlin
      )
    })
    # (3) Apply FDR adjustments across all 186 metabolites
    res_df <- res_df %>%
      mutate(
        EM_beta_FDR = p.adjust(EM_beta_P, method = "BH"),
        EM_FDR_overall = p.adjust(EM_P_overall, method = "BH"),
        EM_FDR_nonlinear = p.adjust(EM_P_nonlinear, method = "BH")
      ) %>%
      # Format P-values as requested
      mutate(across(c(EM_beta_P, EM_beta_FDR, EM_P_overall, EM_P_nonlinear, EM_FDR_overall, EM_FDR_nonlinear),
                    ~ ifelse(.x < 0.001, "<0.001", sprintf("%.3f", .x)),
                    .names = "{.col}_format"))
    return(res_df)
  })
  # Name the list elements by exposure variables
  names(results_list) <- exposures
  return(results_list)
}

# 2.Execute Main & Sensitivity Analyses####

# Main Analysis
Adip_all <- run_EM_analysis(data = metsplc, exposures = exp_vars, outcomes = met_vars, covs = covariates)

# Sens 1: Controls only
Adip_all_ctrl <- run_EM_analysis(data = subset(metsplc, status == 0), exposures = exp_vars, outcomes = met_vars, covs = covariates)

# Sens 2: Person-years >= 2
Adip_all_2y <- run_EM_analysis(data = subset(metsplc, personyear >= 2), exposures = exp_vars, outcomes = met_vars, covs = covariates)

# Sens 3: Excluding Chronic Liver Disease (HLD == 0)
Adip_all_HLD <- run_EM_analysis(data = subset(metsplc, HLD == 0), exposures = exp_vars, outcomes = met_vars, covs = covariates)

