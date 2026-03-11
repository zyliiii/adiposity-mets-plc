# ==============================================================================
# Script Name: 04_Intersection_and_Sensitivity_Filtering.R
# Purpose:     1. Merge E-M and M-D association results
#              2. Filter significant E-M and M-D associations
#              3. Apply strict sensitivity analysis filters
#              4. Extract the final robust intersecting metabolites
# ==============================================================================

# Load required packages
library(dplyr)
library(purrr)

# Define exposure variables
exp_vars <- c("zBMI", "zWaist", "zWHR", "zWHtR", "zABSI", "zHip", "zweightdiff")

# 1. Base Merging & Initial Significance Filtering####

# Identify significant M-D associations overall
CLR_sig <- CLR_all %>% 
  filter(MD_contOR_FDR < 0.05 | MD_P_nonlinear < 0.05)

# Generate lists for merged and significant results
Merge_all <- list()
Adip_sig <- list()
Merge_sig <- list()

for (expo in exp_vars) {
  # Merge all E-M results with M-D results by metabolite
  merged_df <- merge(Adip_all[[expo]], CLR_all, by = "metabolite")
  Merge_all[[expo]] <- merged_df
  # Filter for significant E-M associations
  sig_EM <- merged_df %>%
    filter(EM_beta_FDR < 0.05 | EM_P_nonlinear < 0.05)
  Adip_sig[[expo]] <- sig_EM
  # Filter for significant M-D associations
  Merge_sig[[expo]] <- sig_EM %>%
    filter(MD_contOR_FDR < 0.05 | MD_P_nonlinear < 0.05)
}

# 2. Sensitivity Analysis Filtering on E-M Associations####
# Filter Adip_sig strictly based on the 3 sensitivity analyses

Adip_retained <- purrr::map(exp_vars, function(expo) {
  # Extract main and sensitivity results
  main_df <- Adip_sig[[expo]]
  ctrl_df <- Adip_all_ctrl[[expo]]
  twoy_df <- Adip_all_2y[[expo]]
  HLD_df  <- Adip_all_HLD[[expo]]
  # Join essential columns from sensitivity analyses
  joined_df <- main_df %>%
    left_join(select(ctrl_df, metabolite, ctrl_EM_beta = EM_beta, ctrl_EM_P_nonlinear = EM_P_nonlinear), by = "metabolite") %>%
    left_join(select(twoy_df, metabolite, twoy_EM_beta = EM_beta, twoy_EM_P_nonlinear = EM_P_nonlinear), by = "metabolite") %>%
    left_join(select(HLD_df,  metabolite, HLD_EM_beta = EM_beta,  HLD_EM_P_nonlinear = EM_P_nonlinear),  by = "metabolite")
  # Apply filtering logic:
  # 1. If linear (EM_P_nonlinear >= 0.05): Beta signs must match across all 3 sensitivity tests.
  # 2. If non-linear (EM_P_nonlinear < 0.05): Non-linear P-values must be < 0.05 across all 3 sensitivity tests.
  retained_df <- joined_df %>%
    filter(
      (EM_P_nonlinear >= 0.05 & 
         sign(EM_beta) == sign(ctrl_EM_beta) & 
         sign(EM_beta) == sign(twoy_EM_beta) & 
         sign(EM_beta) == sign(HLD_EM_beta)) |
        (EM_P_nonlinear < 0.05 & 
           ctrl_EM_P_nonlinear < 0.05 & 
           twoy_EM_P_nonlinear < 0.05 & 
           HLD_EM_P_nonlinear < 0.05)
    ) %>%
    # Remove the temporary sensitivity columns to keep the dataset clean
    select(-starts_with("ctrl_"), -starts_with("twoy_"), -starts_with("HLD_"))
  return(retained_df)
})

# Assign names to the list
names(Adip_retained) <- exp_vars

# 3. Final Intersection (Significant in M-D + Robust in E-M Sensitivity analyses)####
# Filter the retained E-M significant metabolites for M-D significance
Merge_sig_sens <- purrr::map(Adip_retained, function(df) {
  df %>% filter(MD_contOR_FDR < 0.05 | MD_P_nonlinear < 0.05)
})
names(Merge_sig_sens) <- exp_vars

# 4. Consolidate Final Intersection Results into a Single Dataframe####
Intersection <- bind_rows(Merge_sig_sens, .id = "Exposure") %>%
  mutate(id = row_number()) %>%
  select(id, Exposure, everything()) 

