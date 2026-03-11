# ===============================================================================
# Script Name: 00_Data_Preparation.R
# Purpose:     Generate a mock dataset (metsplc) for the NCC study to demonstrate 
#              downstream statistical pipelines (E-M-D relationships).
# Output:      metsplc (dataframe containing 644 observations)
# ===============================================================================

library(dplyr)
library(purrr)

# Set seed for reproducibility
set.seed(2026)

# Define sample size (322 cases + 322 controls)
n_pairs <- 322
n_total <- n_pairs * 2

# 1. Generate base structure####
# (Matched Sets and Outcome)
# Set: factor 1 to 322
# status: 1 = PLC case, 0 = matched control
metsplc_base <- data.frame(
  Set = as.factor(rep(1:n_pairs, each = 2)),
  status = rep(c(1, 0), times = n_pairs)
)

# 2. Generate Covariates####
# Match age closely within sets (± 2 years)
set_age <- rnorm(n_pairs, mean = 55, sd = 8)
mock_startage <- rep(set_age, each = 2) + runif(n_total, min = -1.9, max = 1.9)

metsplc_cov <- metsplc_base %>%
  mutate(
    startage = mock_startage,
    edu = factor(sample(1:3, n_total, replace = TRUE), labels = c("Low", "Medium", "High")),
    inc = factor(sample(1:3, n_total, replace = TRUE), labels = c("Low", "Medium", "High")),
    DIET = rnorm(n_total, mean = 50, sd = 10),
    eversmoker = factor(sample(c(0, 1), n_total, replace = TRUE), labels = c("Never", "Ever")),
    everdrinker = factor(sample(c(0, 1), n_total, replace = TRUE), labels = c("Never", "Ever")),
    MET_total_wk = rnorm(n_total, mean = 1500, sd = 500),
    HLD = sample(c(0, 1), n_total, replace = TRUE, prob = c(0.8, 0.2)),
    Hgallstone = sample(c(0, 1), n_total, replace = TRUE, prob = c(0.9, 0.1)),
    HT2DM = sample(c(0, 1), n_total, replace = TRUE, prob = c(0.85, 0.15)),
    personyear = ifelse(status == 1, runif(n_total, 1, 8), runif(n_total, 5, 15))
  )

# 3. Generate Exposures####
# 7 Anthropometric traits with realistic disease associations
metsplc_exp <- metsplc_cov %>%
  mutate(
    # (1) Non-linear (U-shaped) association for BMI:
    BMI = ifelse(status == 0, 
                 rnorm(n(), mean = 23, sd = 2.5), 
                 ifelse(runif(n()) > 0.3, rnorm(n(), mean = 28, sd = 3), rnorm(n(), mean = 19, sd = 2))),
    # (2) Strong positive associations:
    Waist      = ifelse(status == 0, rnorm(n(), mean = 80, sd = 8), rnorm(n(), mean = 90, sd = 8)),
    WHtR       = ifelse(status == 0, rnorm(n(), mean = 0.48, sd = 0.04), rnorm(n(), mean = 0.54, sd = 0.04)),
    weightdiff = ifelse(status == 0, rnorm(n(), mean = 1, sd = 4), rnorm(n(), mean = 5, sd = 5)),
    # (3) Weak positive associations:
    WHR  = ifelse(status == 0, rnorm(n(), mean = 0.84, sd = 0.07), rnorm(n(), mean = 0.86, sd = 0.07)),
    Hip  = ifelse(status == 0, rnorm(n(), mean = 96, sd = 7), rnorm(n(), mean = 98, sd = 7)),
    ABSI = ifelse(status == 0, rnorm(n(), mean = 0.080, sd = 0.005), rnorm(n(), mean = 0.082, sd = 0.005))
  ) 

# 4. Generate Mediators####
# (186 plasma metabolites with implanted biological signals)
# Create a composite adiposity score to drive the exposure correlation
adip_score <- as.numeric(scale(
  metsplc_exp$BMI + metsplc_exp$Waist + metsplc_exp$WHR + 
    metsplc_exp$WHtR + metsplc_exp$ABSI + metsplc_exp$Hip + metsplc_exp$weightdiff
))


# Initialize an empty matrix for the 186 metabolites
mat_met <- matrix(NA, nrow = n_total, ncol = 186)

for (j in 1:186) {
  if (j <= 20) {
    # Implanted Signal (met1 to met20): 
    # Log-scale is positively driven by Adiposity (Exposure) + Status (Outcome) + Noise
    log_val <- rnorm(n_total, mean = 4, sd = 0.5) + 
      0.5 * adip_score + 
      0.8 * metsplc_exp$status
    # Convert back to natural scale 
    mat_met[, j] <- exp(log_val)
  } else {
    # Random Noise (met21 to met186): 
    # Standard log-normal distribution, no relation to exposure or outcome
    mat_met[, j] <- rlnorm(n_total, meanlog = 4.5, sdlog = 1)
  }
}

# Assign column names (met1 to met186) and bind to the main dataset
colnames(mat_met) <- paste0("met", 1:186)

metsplc_raw <- cbind(metsplc_exp, mat_met)

# ===============================================================================
# 5. Data Transformation (log2 and z-scale for metabolites, z-scale for exposures)####
metsplc <- metsplc_raw %>%
  # (1) Log2 transformation and Z-scaling for all metabolites
  mutate(
    across(
      .cols = starts_with("met"),
      .fns = ~ as.numeric(scale(log2(.x))),
      .names = "zlog2{.col}"
    )
  ) %>%
  # (2) Z-scaling for the 7 exposure variables
  mutate(
    across(
      .cols = c("BMI", "Waist", "WHR", "WHtR", "ABSI", "Hip", "weightdiff"),
      .fns = ~ as.numeric(scale(.x)),
      .names = "z{.col}"
    )
  )


# Clean up environment
rm(metsplc_base, metsplc_cov, metsplc_exp, metsplc_raw, mat_met, set_age, mock_startage, n_pairs, n_total, adip_score, j, log_val)
