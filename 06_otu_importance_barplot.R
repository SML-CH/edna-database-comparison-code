# ============================================================
# Visualize OTU importance scores with taxonomic annotation
# For CMBL WQI model
# ============================================================

# Load required packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# Step 1: Copy the table from PDF (including header)
# Select the table content in PDF, press Ctrl+C (Windows) or Command+C (Mac)

# Step 2: Read data from clipboard
otu_table <- read.delim("clipboard", header = TRUE, sep = "\t", row.names = NULL)

# For Mac users, if the above fails, use this line instead (comment out the line above)
# otu_table <- read.delim(pipe("pbpaste"), header = TRUE)

# Check the first few rows of data
head(otu_table)

# Step 3: Extract required columns, use only Species for Latin names
df <- otu_table %>%
  select(OTU, Mean.importance.score, Class, Species) %>%
  rename(Importance = Mean.importance.score) %>%
  # Use Species directly as Latin name
  mutate(Latin_name = as.character(Species))

# Step 4: Sort by importance
df <- df %>%
  arrange(Importance) %>%
  mutate(OTU = factor(OTU, levels = OTU))

# Step 5: High-impact journal color scheme (based on Class)
class_present <- unique(df$Class)

# Nature-style colors for different classes
nature_colors <- c(
  "Gastropoda" = "#E64B35",      # Gastropoda - Coral red
  "Bivalvia" = "#4DBBD5",        # Bivalvia - Cyan blue
  "Insecta" = "#00A087",         # Insecta - Dark green
  "Clitellata" = "#3C5488",      # Clitellata - Deep blue
  "Malacostraca" = "#F39B7F",    # Malacostraca - Orange
  "norank" = "#8491B4",          # Unclassified - Gray blue
  "Neogastropoda" = "#DC0000",   # Neogastropoda - Red
  "Unionida" = "#7E6148",        # Unionida - Brown
  "Veneroida" = "#B09C85"        # Veneroida - Beige
)

# Keep only classes present in the data
class_colors <- nature_colors[names(nature_colors) %in% class_present]

# Automatically assign colors for new classes not in the preset list
new_class <- class_present[!class_present %in% names(nature_colors)]
if (length(new_class) > 0) {
  extra_colors <- setNames(
    brewer.pal(min(length(new_class), 8), "Set2"),
    new_class
  )
  class_colors <- c(class_colors, extra_colors)
}

# Step 6: Plot (Latin names only show Species, placed to the right of bars)
p <- ggplot(df, aes(x = Importance, y = OTU, fill = Class)) +
  geom_col(width = 0.8) +
  geom_text(aes(x = Importance, label = Latin_name), 
            hjust = -0.05,
            size = 3,
            fontface = "italic",
            family = "serif") +
  scale_fill_manual(values = class_colors) +
  labs(x = "Mean Importance Score", 
       y = "", 
       title = "",
       fill = "Class") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8, face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_line(color = "gray90"),
    plot.margin = margin(t = 5, r = 70, b = 5, l = 5, unit = "pt")
  ) +
  scale_x_continuous(
    limits = c(0, max(df$Importance) * 1.45),
    expand = c(0, 0)
  )

# Display the plot
print(p)

# Automatically set compact plot height based on number of OTUs
plot_height <- max(6, nrow(df) * 0.22)

# Check current working directory
getwd()

# Save high-resolution images
ggsave("OTU_importance_CMBL_WQI.png", 
       plot = p,
       width = 12, height = plot_height, dpi = 300)

ggsave("OTU_importance_CMBL_WQI.pdf", 
       plot = p,
       width = 12, height = plot_height, dpi = 300)


####  boxplot  ####
#### Single group boxplot ####
# Import data
data <- read.delim("clipboard", header = TRUE)

# Extract value column
x <- data$value

# Draw vertical boxplot
boxplot(x,
        horizontal = FALSE,   # Vertical
        col = "white",
        border = "black",
        main = "",
        ylab = "")

# Calculate mean
mean_x <- mean(x, na.rm = TRUE)

# Add mean as red plus sign (bold)
points(1, mean_x, col = "red", pch = 3, cex = 1.8, lwd = 2)

# Label mean value next to it (red, right-aligned)
text(1.2, mean_x, round(mean_x, 2), col = "red")

# Calculate five-number summary (min, Q1, median, Q3, max)
qs <- fivenum(x)

# Label five-number summary on the left
text(0.8, qs, labels = round(qs, 2), adj = 1)

#### Two-group boxplot ####
# Set working directory
setwd("E:/PhD_Year1-2/Database_comparison/eDNA_submission/Revision_results")

# Import data
data <- read.csv("Species_overlap_rate_boxplot.csv", header = TRUE)

# Boxplot + points + high-impact journal style
p <- ggplot(data, aes(x = Group, y = Value)) +
  # Boxplot, black outline, no fill, bold whisker lines
  geom_boxplot(width = 0.6, fill = NA, color = "black", outlier.shape = NA, size = 1) +
  # Jittered points
  geom_jitter(aes(color = Group), width = 0.15, size = 2, alpha = 0.8) +
  # Point colors
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  # Theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(size = 1)  # Bold axes
  ) +
  # Significance annotation
  stat_compare_means(method = "t.test", label = "p.signif", size = 5) +
  labs(x = "", y = "") +
  # Add full border
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1.2))

# Display the plot
print(p)
