# eDNA Database Comparison Code

R code accompanying the eDNA database comparison manuscript. The repository contains ten analyses covering community similarity, family-level overlap, random-forest prediction of biomonitoring indices, regression, rarefaction, eDNA–morphology congruence, watershed-level validation, and GLM analysis of watershed and year effects.

## Analyses

1. **`01_bray_curtis_similarity.R`** — Site-level Bray–Curtis similarity between morphological and eDNA family-level communities.
2. **`02_family_venn_diagram.R`** — Area-proportional family-level Euler/Venn comparison of CMBL, NCBI, and morphology datasets.
3. **`03_css_loocv_random_forest.R`** — Leakage-free CSS-normalized LOOCV random forest and feature-number sensitivity analysis.
4. **`04_regression.R`** — Regression analysis of observed versus predicted biomonitoring index values.
5. **`05_site_rarefaction_curves.R`** — Site-level rarefaction curves for unannotated, CMBL, and NCBI OTU matrices.
6. **`06_otu_importance_barplot.R`** — Important-OTU visualization and species-overlap plotting.
7. **`07_eDNA–Morphology Congruence across Databases.R`** — Comparison of low- and high-congruence sites across CMBL and NCBI databases.
8. **`08_lowocv_random_forest.R`** — Leave-one-watershed-out random-forest validation with training-set CSS normalization and feature selection.
9. **`09_GLM Analysis of Watershed and Year Effects`** — GLM analysis of watershed and year effects on observed and predicted WQI/BMWP after excluding Qiantang River samples.
10. **`10_rarefaction_loocv_random_forest`** — Leakage-free rarefaction-normalized LOOCV random forest for BMWP prediction.

## Data and paths

Raw data are not included. Before running a code file:

1. inspect the input file names declared near the beginning of the script;
2. place the required data in a local analysis folder;
3. replace the example `setwd()` path with your local path;
4. keep generated results separate from the original input data.

Code 9 prints its results to the R console and intentionally does not generate output files.

## Installation

Common CRAN dependencies can be installed with:

```r
install.packages(c(
  "caret", "cowplot", "dplyr", "e1071", "eulerr", "ggplot2",
  "ggpubr", "ggrepel", "gridExtra", "patchwork", "purrr",
  "RColorBrewer", "ranger", "readr", "readxl", "stringr",
  "tidyr", "tidyverse", "vegan", "viridis"
))
```

CSS-normalization scripts additionally require:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Biobase", "metagenomeSeq"))
```

## Running a code file

From R:

```r
source("01_bray_curtis_similarity.R")
```

From a terminal:

```sh
Rscript 01_bray_curtis_similarity.R
```

See **`R CODE STATEMENT`** for detailed descriptions and reproducibility notes for Codes 1–10.

