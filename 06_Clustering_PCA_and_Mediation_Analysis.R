# ==============================================================================
# Script Name: 06_Clustering_PCA_and_Mediation_Analysis.R
# Purpose:     1. Hierarchical clustering of significant metabolites
#              2. Extract PC1 for each cluster
#              3. Mediation analysis using Hayes PROCESS macro
# ==============================================================================

# Load required packages
library(dplyr)
library(ClustOfVar)

# 1. Data Preparation for PROCESS Macro####
PROCESS_data <- metsplc

# Convert factor variables to strict numeric (0/1)
PROCESS_data <- PROCESS_data %>%
  mutate(
    status = as.numeric(as.character(status)),
    HLD = as.numeric(as.character(HLD)),
    Hgallstone = as.numeric(as.character(Hgallstone)),
    HT2DM = as.numeric(as.character(HT2DM)),
    # Convert 'ever/never' factors to 1/0
    eversmoker = ifelse(eversmoker == "Ever", 1, 0),
    everdrinker = ifelse(everdrinker == "Ever", 1, 0),
    # Create Dummy Variables for Categorical Covariates
    edu_2 = ifelse(edu == "Medium", 1, 0),
    edu_3 = ifelse(edu == "High", 1, 0),
    inc_2 = ifelse(inc == "Medium", 1, 0),
    inc_3 = ifelse(inc == "High", 1, 0)
  )

# 2. Hierarchical Clustering and PC1 Extraction####

# Pre-defined optimal number of clusters (k) for each exposure (based on bootstrapped-meanCR)
k_map <- list(
  "zBMI" = 2, 
  "zWaist" = 5, 
  "zWHR" = 2, 
  "zWHtR" = 3, 
  "zABSI" = 6, 
  "zHip" = 4, 
  "zweightdiff" = 5
)

# List to store all new PC1 variables
pca_variables_list <- list()
# List to store clustering models
Clusters_res <- list()

for (expo in names(k_map)) {
  k <- k_map[[expo]]
  # Extract significant metabolites for the current exposure
  zlog2_vars <- Intersection %>% 
    filter(Exposure == expo) %>% 
    pull(metabolite)
  # Subset data for clustering
  var_data <- PROCESS_data[, zlog2_vars, drop = FALSE]
  # Perform Hierarchical Clustering
  tree <- hclustvar(X.quanti = var_data)
  cluster_cut <- cutreevar(tree, k = k)
  # Save cluster grouping info
  Clusters_res[[expo]] <- cluster_cut
  # Extract PCA scores (PC1 for each cluster)
  pca_scores <- as.data.frame(cluster_cut$scores)
  # Rename columns to standard format
  colnames(pca_scores) <- paste0(expo, "_PC1_", 1:k)
  # Bind all generated PC1 columns back to the main PROCESS_data
  PROCESS_data <- cbind(PROCESS_data, pca_scores)
}


# 3. Mediation Analysis (Hayes PROCESS Macro)####
# NOTE: Make sure the PROCESS macro (version 4.3) is loaded before running!!
# The PROCESS.R can be download at processmacro.org.

# Define covariates
PROCESS_covs <- c("startage", "edu_2", "edu_3", "inc_2", "inc_3", 
                  "eversmoker", "everdrinker", "DIET", "MET_total_wk", 
                  "HLD", "Hgallstone", "HT2DM")

# List to store PROCESS output text
PROCESS_results <- list()

for (expo in names(k_map)) {
  k <- k_map[[expo]]
  # Define X and M variable names
  x_var <- paste0(expo)
  m_vars <- paste0(x_var, "_PC1_", 1:k)
  
  cat("Running PROCESS macro for:", expo, "...\n")
  
  # Capture the printed output of the process() function
  output_text <- capture.output({
    process(
      data     = PROCESS_data, 
      y        = "status", 
      x        = x_var, 
      m        = m_vars, 
      cov      = c("startage","edu_2","edu_3","inc_2","inc_3","eversmoker","everdrinker",
                   "DIET","MET_total_wk","HLD","Hgallstone","HT2DM"), 
      model    = 4, 
      boot     = 2000, 
      seed     = 666, 
      normal   = 1, 
      conf     = 95, 
      xmtest   = 1,
      mcx      = 0, 
      total    = 0, 
      contrast = 0, 
      modelres = 0    )
  })
  # Save the result string in the list
  PROCESS_results[[expo]] <- output_text
}


# To view the results for a specific exposure (e.g., BMI), run the following:
cat(PROCESS_results[["zBMI"]], sep = "\n")

