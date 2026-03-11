# ==============================================================================
# Script Name: 01_Baseline_Characteristics.R
# Purpose:     1. Descriptive statistics of adiposity exposures
#              2. Spearman correlation matrix among exposures
#              3. Baseline characteristics comparison by PLC status
# ==============================================================================

# Load required packages
library(dplyr)
library(purrr)
library(writexl)
library(moonBook)

# Define exposure variables and their corresponding display labels
exp_vars <- c("BMI", "Waist", "WHR", "WHtR", "ABSI", "Hip", "weightdiff")
exp_labels <- c("BMI", "WC", "WHR", "WHtR", "ABSI", "HC", "Adult weight gain")

# 1. Descriptive Statistics of Adiposity Exposures####
descriptive_adiposity <- map2_dfr(exp_vars, exp_labels, function(var, label) {
  x <- metsplc[[var]]
  data.frame(
    Var = label,
    min = min(x, na.rm = TRUE),
    x   = mean(x, na.rm = TRUE),
    SD  = sd(x, na.rm = TRUE),
    P25 = quantile(x, 0.25, na.rm = TRUE),
    M   = median(x, na.rm = TRUE),
    P75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
})



# 2. Spearman Correlation Between Adiposity Indicators####
cor_matrix_df <- local({
  n <- length(exp_vars)
  mat <- matrix("", n, n, dimnames = list(exp_vars, exp_vars))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      res <- cor.test(metsplc[[exp_vars[i]]], metsplc[[exp_vars[j]]], method = "spearman")
      mat[j, i] <- sprintf("%.2f%s", res$estimate, ifelse(res$p.value < 0.05, "*", ""))
    }
  }
  as.data.frame(mat)
})



# 3. Baseline Characteristics Comparison by Status####
# Ensure status is treated as a factor
metsplc$status_factor <- as.factor(metsplc$status)

# Generate baseline table
table_baseline <- mytable(
  status_factor ~ startage + edu + inc + eversmoker + everdrinker + 
    HLD + Hgallstone + HT2DM + 
    MET_total_wk + DIET + 
    BMI + Waist + WHR + WHtR + ABSI + Hip + weightdiff,
  data       = metsplc,
  method     = 3, 
  catMethod  = 0, 
  digits     = 2, 
  show.all   = TRUE, 
  show.total = TRUE
)

table_baseline
