# adiposity-met-plc

This repository now contains the R scripts used for the nested case-control (NCC) analysis portion of our study published in PLOS Medicine (https://doi.org/10.1371/journal.pmed.1004910). 


## Overview of the Analysis Pipeline
The analysis pipeline is modularized into 7 sequential R scripts. To comply with data privacy and ethical regulations, the raw clinical dataset is not publicly available. However, `00_Data_Preparation.R` includes code to generate a statistically representative mock dataset, allowing researchers to fully reproduce and test the analytical pipeline.


### Scripts Description:
* **`00_Data_Preparation.R`**: Generates a mock dataset with simulated structural and biological associations. Applies log2-transformation and Z-scaling for metabolic features, and Z-scaling for exposure traits.
* **`01_Baseline_Characteristics.R`**: Summarizes the baseline characteristics of the study population.
* **`02_Exposure_Mediator_Association.R`**: Performs linear and non-linear association analyses between adiposity indicators (Exposures) and plasma metabolites (Mediators), incorporating strict FDR adjustments.
* **`03_Mediator_Disease_Association.R`**: Conducts conditional logistic regression (CLR) to assess the associations between metabolites (Mediators) and the risk of primary liver cancer (Disease).
* **`04_Intersection_and_Sensitivity_Filtering.R`**: Extracts the robustly significant intersection of metabolites that passed strict multi-level sensitivity analyses.
* **`05_Network_Analysis.R`**: Constructs a data-driven causal network skeleton using the PC algorithm and identifies densely connected metabolite modules via the Walktrap algorithm.
* **`06_Clustering_PCA_and_Mediation_Analysis.R`**: Performs hierarchical clustering on significant metabolites, extracts the first principal component (PC1) for each cluster, and conducts mediation analysis utilizing the Hayes PROCESS macro.


## Requirements
* **R** (version >= 4.4.1 recommended)
* **Required Packages:** `dplyr`, `purrr`, `survival`, `pcalg`, `igraph`, `ClustOfVar`.
* **PROCESS Macro:** The mediation analysis requires the Hayes PROCESS macro for R (`process.R`), which should be loaded into the environment prior to running script `06`.


## Citation
If you use these codes in your research, please cite our paper:
Li Z-Y, Li H-L, Wang J, Shen Q-M, Zou Y-X, Yang D-N, et al. Metabolomic insights into associations between adiposity markers and liver cancer risk: Results from a prospective cohort study and Mendelian randomization analysis. PLoS Med. 2026;23: e1004910. doi:10.1371/journal.pmed.1004910


Thanks for your interest! 
The code for the two-sample MR analysis is under preparation and will be available soon.
