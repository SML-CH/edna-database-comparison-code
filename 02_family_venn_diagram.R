############################################################
## CMBL vs NCBI vs Morphology family-level comparison
## Venn / Euler diagram with area-proportional scaling
############################################################

# Load required packages
library(eulerr)
library(dplyr)
library(grid)

# Set working directory
setwd("E:/PhD/Database_comparison")

# =========================
# 1. Read data
# =========================
set1 <- read.csv("CMBL_family_level_data.csv", header = TRUE, stringsAsFactors = FALSE)
set2 <- read.csv("NCBI_family_level_data.csv", header = TRUE, stringsAsFactors = FALSE)
set3 <- read.csv("Zhejiang_morphological_family_level_data.csv", header = TRUE, stringsAsFactors = FALSE)

# =========================
# 2. Extract family names
#    Remove NA, empty values, and duplicates
# =========================
CMBL_families <- set1$family %>%
  na.omit() %>%
  trimws() %>%
  unique()

NCBI_families <- set2$family %>%
  na.omit() %>%
  trimws() %>%
  unique()

Morphology_families <- set3$family %>%
  na.omit() %>%
  trimws() %>%
  unique()

# =========================
# 3. Create list of sets
# =========================
family_data <- list(
  CMBL = CMBL_families,
  NCBI = NCBI_families,
  Morphology = Morphology_families
)

# =========================
# 4. Fit area-proportional Euler / Venn diagram
# =========================
fit <- euler(family_data)

# View fitting performance
print(fit)

# =========================
# 5. High-impact journal style plotting
# =========================
plot(
  fit,
  
  # Fill colors: soft, low saturation, publication-ready
  fills = list(
    fill = c("#E68D3D", "#E26472", "#6270B7"),
    alpha = 0.55
  ),
  
  # Circle borders
  edges = list(
    col = "grey25",
    lwd = 1.2
  ),
  
  # Set labels
  labels = list(
    font = 2,
    fontsize = 15,
    col = "grey10"
  ),
  
  # Region quantities (counts)
  quantities = list(
    fontsize = 14,
    font = 2,
    col = "grey10"
  ),
  
  # Do not show percentages, only counts
  legend = FALSE,
  
  # Background
  main = NULL
)

############################################################
## Export high-resolution figure
############################################################

pdf(
  file = "Family_level_area_proportional_Venn.pdf",
  width = 7,
  height = 6
)

plot(
  fit,
  fills = list(
    fill = c("#E68D3D", "#E26472", "#6270B7"),
    alpha = 0.55
  ),
  edges = list(
    col = "grey25",
    lwd = 1.2
  ),
  labels = list(
    font = 2,
    fontsize = 15,
    col = "grey10"
  ),
  quantities = list(
    fontsize = 14,
    font = 2,
    col = "grey10"
  ),
  legend = FALSE,
  main = NULL
)

dev.off()
