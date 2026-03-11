# ===================================================================================
# Script Name: 03_Mediator_Disease_Association.R
# Purpose:     Conditional Logistic Regression for Mediator-Disease (M-D) association
#              Includes Main analysis and Sensitivity analysis
# ===================================================================================

# Load required packages
library(survival)
library(dplyr)
library(purrr)

# Define exposure variables
met_vars <- paste0("zlog2met", 1:186)
# Define covariates
MD_covars <- "startage + eversmoker + everdrinker + BMI + MET_total_wk + DIET + HLD + Hgallstone + HT2DM"

# 1. Main Analysis: Conditional Logistic Regression####
run_MD_clogit <- function(data, outcomes, covs) {
  res_df <- purrr::map_dfr(outcomes, function(outc) {
    # Linear clogit
    form_lin <- as.formula(paste("status ~", outc, "+", covs, "+ strata(Set)"))
    mod_lin <- clogit(form_lin, data = data, method = "exact")
    stats_lin <- summary(mod_lin)$coefficients
    
    coef_val <- stats_lin[outc, "coef"]
    se_val <- stats_lin[outc, "se(coef)"]
    p_lin <- stats_lin[outc, "Pr(>|z|)"]
    
    or_lin <- exp(coef_val)
    lci_lin <- exp(coef_val - 1.96 * se_val)
    uci_lin <- exp(coef_val + 1.96 * se_val)
    # RCS clogit
    form_rcs <- as.formula(paste("status ~ rcs(", outc, ", 3) +", covs, "+ strata(Set)"))
    mod_rcs <- clogit(form_rcs, data = data)
    stats_rcs <- summary(mod_rcs)$coefficients
    # Exact matching of coefficient names to extract P-values
    p_overall <- stats_rcs[paste0("rcs(", outc, ", 3)", outc), "Pr(>|z|)"]
    p_nonlin <- stats_rcs[paste0("rcs(", outc, ", 3)", outc, "'"), "Pr(>|z|)"]
    
    data.frame(
      metabolite = outc,
      MD_contOR = or_lin,
      MD_contOR_LCI = lci_lin,
      MD_contOR_UCI = uci_lin,
      MD_contOR_ci = sprintf("%.2f (%.2f, %.2f)", or_lin, lci_lin, uci_lin),
      MD_contOR_P = p_lin,
      MD_P_overall = p_overall,
      MD_P_nonlinear = p_nonlin
    )
  })
  # Apply FDR and formatting
  res_df <- res_df %>%
    mutate(
      MD_contOR_FDR = p.adjust(MD_contOR_P, method = "BH"),
      MD_FDR_overall = p.adjust(MD_P_overall, method = "BH"),
      MD_FDR_nonlinear = p.adjust(MD_P_nonlinear, method = "BH")
    ) %>%
    mutate(across(c(MD_contOR_P, MD_contOR_FDR, MD_P_overall, MD_P_nonlinear, MD_FDR_overall, MD_FDR_nonlinear),
                  ~ ifelse(.x < 0.001, "<0.001", sprintf("%.3f", .x)),
                  .names = "{.col}_format"))
  return(res_df)
}

# Execute Main Analysis
CLR_all <- run_MD_clogit(data = metsplc, outcomes = met_vars, covs = MD_covars)

# 2. Sensitivity Analysis: Unconditional Logistic Regression (Subset without HLD)####
# Define covariates
MD_covars_sens <- "startage + eversmoker + everdrinker + BMI + MET_total_wk + DIET + Hgallstone + HT2DM" 

run_MD_glm <- function(data, outcomes, covs) {
  res_df <- purrr::map_dfr(outcomes, function(outc) {
    form_glm <- as.formula(paste("status ~", outc, "+", covs))
    mod_glm <- glm(form_glm, data = data, family = binomial)
    
    coef_est <- coef(summary(mod_glm))

    ci <- suppressMessages(confint(mod_glm, level = 0.95)[outc, ])
    
    or_val <- exp(coef_est[outc, "Estimate"])
    ci_lower <- exp(ci[1])
    ci_upper <- exp(ci[2])
    p_val <- coef_est[outc, "Pr(>|z|)"]
    
    data.frame(
      metabolite = outc,
      OR = or_val,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      OR_CI = sprintf("%.2f (%.2f, %.2f)", or_val, ci_lower, ci_upper),
      P_value = p_val
    )
  })
  # Apply FDR and formatting
  res_df <- res_df %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    mutate(across(c(P_value, FDR),
                  ~ ifelse(.x < 0.001, "<0.001", sprintf("%.3f", .x)),
                  .names = "{.col}_format"))
  return(res_df)
}

# Execute Sensitivity Analysis
Sens_GLM_all <- run_MD_glm(data = subset(metsplc, HLD == 0), outcomes = met_vars, covs = MD_covars_sens)

