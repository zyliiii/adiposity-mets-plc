# ==============================================================================
# Script Name: 05_Network_Analysis.R
# Purpose:     1. Construct a data-driven network using the PC algorithm
#              2. Identify densely connected metabolite modules via Walktrap
# ==============================================================================

# Load required packages
library(pcalg)
library(igraph)
library(dplyr)

# 1. Data Preparation####

# Extract unique significant metabolites from the Intersection dataset
unique_mets <- Intersection %>%
  distinct(metabolite, .keep_all = TRUE) %>%
  pull(metabolite)

# Subset the main dataset (metsplc) using these unique metabolites
var_data <- metsplc[, unique_mets]
var_data <- var_data[,1:15]

# Calculate correlation matrix and sample size required for the PC algorithm
corMatrix <- cor(var_data, use = "pairwise.complete.obs")
n_samples <- nrow(var_data)

# Create the sufficient statistics list required by pcalg::pc()
suffStat <- list(C = corMatrix, n = n_samples)

# 2. PC Algorithm (Causal Network Skeleton Construction)####

# Apply the PC algorithm using Gaussian Conditional Independence test
pc_result <- pc(
  suffStat = suffStat, 
  indepTest = gaussCItest, 
  alpha = 0.05, 
  labels = colnames(var_data)
)

# Print PC network summary
print(pc_result)

# 3. Walktrap Algorithm (Module Identification)####

# Extract adjacency matrix from the PC result and convert to an igraph object
adjacency_matrix <- as(pc_result@graph, "matrix")
net_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "max")

# Run Walktrap clustering
wt_clusters <- cluster_walktrap(net_graph, steps = 4)

# 4. Results Extraction & Visualization####

# Extract cluster membership for each metabolite
membership_vec <- membership(wt_clusters)

# Create a dataframe linking each metabolite node to its identified cluster
cluster_res <- data.frame(
  node = names(membership_vec),
  cluster = as.vector(membership_vec),
  stringsAsFactors = FALSE
)


# Plot the network graph with identified clusters
set.seed(666)
plot(wt_clusters, net_graph, 
     vertex.size = 15, 
     vertex.label.cex = 0.8)

