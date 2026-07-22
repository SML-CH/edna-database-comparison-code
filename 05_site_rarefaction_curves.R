############################################################
## OTU rarefaction curves for three data matrices
## Raw OTUs, CMBL-annotated OTUs, and NCBI-annotated OTUs
############################################################

# Load required package
library(vegan)

# Set working directory
setwd("E:/PhD/Database_comparison/eDNA_revision_results")

# Function to read OTU table and format for rarefaction curves
read_otu_table <- function(file) {
  otu <- read.csv(file, row.names = 1, check.names = FALSE)
  otu <- as.matrix(otu)
  storage.mode(otu) <- "numeric"
  otu[is.na(otu)] <- 0
  otu <- round(otu)
  
  # rarecurve requires rows = samples, columns = OTUs
  otu_site <- t(otu)
  
  # Remove empty samples and empty OTUs
  otu_site <- otu_site[rowSums(otu_site) > 0, , drop = FALSE]
  otu_site <- otu_site[, colSums(otu_site) > 0, drop = FALSE]
  
  return(otu_site)
}

# Function to plot rarefaction curve
plot_rarecurve <- function(otu_site, title_text, line_col) {
  rarecurve(
    otu_site,
    step = 1000,
    label = FALSE,
    col = adjustcolor(line_col, alpha.f = 0.35),
    xlab = "Sequencing depth (reads)",
    ylab = "Observed OTUs",
    main = title_text,
    lwd = 0.8
  )
}

# Read data
otu_raw  <- read_otu_table("rarefaction_raw.csv")
otu_cmbl <- read_otu_table("rarefaction_cmbl.csv")
otu_ncbi <- read_otu_table("rarefaction_ncbi.csv")

# Export to PDF
pdf("Figure_Sx_OTU_rarefaction_curves_three_matrices.pdf", width = 15, height = 5)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

plot_rarecurve(otu_raw,  "Taxonomy-free raw OTUs", "#3B82F6")
plot_rarecurve(otu_cmbl, "CMBL-annotated OTUs", "#10B981")
plot_rarecurve(otu_ncbi, "NCBI-annotated OTUs", "#F59E0B")

dev.off()
