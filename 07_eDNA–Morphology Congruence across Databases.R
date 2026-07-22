############################################################
## Comparison of low vs high congruence between eDNA and morphology
## For both CMBL and NCBI databases
############################################################

# Load required packages
library(ggplot2)
library(ggpubr)
library(patchwork)

# Set directory
setwd("E:/PhD/Database_comparison/eDNA_revision_results/Congruence_analysis")

# Read data
CMBL <- read.csv("CMBL_congruence_analysis.csv", check.names = FALSE)
NCBI <- read.csv("NCBI_congruence_analysis.csv", check.names = FALSE)

# Prepare data: group names and ordering
prepare_data <- function(data) {
  data$group <- factor(
    data$group,
    levels = c("low-congruence", "high-congruence"),
    labels = c("Low", "High")
  )
  data$Reads_log10 <- log10(data$Reads)
  data
}

CMBL <- prepare_data(CMBL)
NCBI <- prepare_data(NCBI)

# Color scheme
group_colors <- c(
  "Low"  = "#4C78A8",
  "High" = "#E45756"
)

# Function to create individual panel
make_panel <- function(data, variable, title) {
  
  ggplot(data, aes(x = group, y = .data[[variable]], fill = group)) +
    geom_boxplot(
      width = 0.55,
      linewidth = 0.75,
      outlier.shape = NA,
      alpha = 0.90
    ) +
    geom_jitter(
      aes(color = group),
      width = 0.13,
      size = 2,
      alpha = 0.65
    ) +
    stat_compare_means(
      comparisons = list(c("Low", "High")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.02,
      bracket.size = 0.6,
      size = 5
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_y_continuous(
      expand = expansion(mult = c(0.04, 0.22))
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_classic(base_size = 15) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        size = 15,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 12,
        color = "black"
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.7
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = margin(8, 5, 8, 5)
      ),
      plot.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.45
      ),
      plot.margin = margin(6, 6, 6, 6)
    )
}

# First row: CMBL
p_a1 <- make_panel(
  CMBL, "Reads_log10",
  "Sequencing reads\nCount (log10)"
)

p_a2 <- make_panel(
  CMBL, "OTUs",
  "eDNA\nOTU Richness"
)

p_a3 <- make_panel(
  CMBL, "eDNA Total Families",
  "eDNA\nFamily Richness"
)

p_a4 <- make_panel(
  CMBL, "Morphology Species Count",
  "Morphology\nOTU Richness"
)

p_a5 <- make_panel(
  CMBL, "Morphology Family Count",
  "Morphology\nFamily Richness"
)

# Second row: NCBI
p_b1 <- make_panel(
  NCBI, "Reads_log10",
  "Sequencing reads\nCount (log10)"
)

p_b2 <- make_panel(
  NCBI, "OTUs",
  "eDNA\nOTU Richness"
)

p_b3 <- make_panel(
  NCBI, "eDNA Total Families",
  "eDNA\nFamily Richness"
)

p_b4 <- make_panel(
  NCBI, "Morphology Species Count",
  "Morphology\nOTU Richness"
)

p_b5 <- make_panel(
  NCBI, "Morphology Family Count",
  "Morphology\nFamily Richness"
)

# Combine plots: a = CMBL, b = NCBI
figure_final <-
  (p_a1 + p_a2 + p_a3 + p_a4 + p_a5) /
  (p_b1 + p_b2 + p_b3 + p_b4 + p_b5) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        size = 25,
        face = "bold"
      )
    )
  ) &
  theme(
    plot.tag.position = c(-0.06, 1.04)
  )

print(figure_final)

# Save as PDF
ggsave(
  filename = file.path(
    root_dir,
    "eDNA_and_morphology_congruence_comparison_CMBL_vs_NCBI.pdf"
  ),
  plot = figure_final,
  width = 13.5,
  height = 9.2,
  units = "in",
  device = cairo_pdf
)

# Save as PNG
ggsave(
  filename = file.path(
    root_dir,
    "eDNA_and_morphology_congruence_comparison_CMBL_vs_NCBI.png"
  ),
  plot = figure_final,
  width = 13.5,
  height = 9.2,
  units = "in",
  dpi = 600,
  bg = "white"
)
