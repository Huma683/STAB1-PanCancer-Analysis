# Pan-Cancer Analysis of STAB1 Expression, Prognosis, and Immune Infiltration

This repository contains R, Bash, and HPC cluster submission scripts used for my M.Sc. Bioinformatics thesis at National Research Tomsk State University. The pipeline evaluates *STAB1* mRNA expression, patient survival outcomes, and tumor-associated macrophage (TAM) infiltration across 33 TCGA cancer types, with external validation using GEO and Human Protein Atlas (HPA) datasets.

## Repository Structure

### Survival & Expression Analysis (R)
* `stab1_baseR_survival.R` — Evaluates overall survival (OS) and disease-free survival (DFS) across TCGA cohorts.
* `stab1_baseR_full_survival_stage.R` — Performs stage-stratified survival modeling and Cox regression analyses.
* `stab1_visualization_baseR.R` — Generates pan-cancer expression boxplots, forest plots, and Kaplan-Meier survival curves.

### External & Immune Validation (R)
* `stab1_geo_validation_baseR.R` — Validates transcriptomic findings across independent GEO cohorts.
* `stab1_hpa_validation_baseR.R` — Cross-references mRNA trends with Human Protein Atlas (HPA) protein-level evidence.
* `stab1_immune_validation_baseR.R` — Assesses correlations between *STAB1* expression and immune cell infiltration.

### Automation & HPC Workflows (Bash / Cluster)
* `run_stab1_full_baseR.sh` & `stab1_full_baseR.sub` — Automates batch execution on a Linux HPC cluster.
* `run_stab1_full_transfer.sh` & `stab1_full_transfer.sub` — Manages automated data transfer and result aggregation.
