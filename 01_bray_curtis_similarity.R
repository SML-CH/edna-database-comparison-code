##############################################
#  Bray-Curtis Similarity Analysis (using union of all families)
##############################################

# Load required packages
library(vegan)
library(tidyverse)

setwd("E:/PhD_Year/Database_comparison/eDNA_submission/Revision_results")
# Read data
morph_data <- read.csv("Zhejiang_morphological_family_level_data.csv", row.names = 1, stringsAsFactors = FALSE)
edna_data  <- read.csv("NCBI_family_level_data.csv", row.names = 1, stringsAsFactors = FALSE)

cat("Morphological data dimensions: ", paste(dim(morph_data), collapse = " x "), "\n")
cat("eDNA data dimensions: ", paste(dim(edna_data), collapse = " x "), "\n")

##############################################
##############################################
# Use the union of all families
##############################################
cat("\n=== Constructing union of families ===\n")

all_taxa <- union(rownames(morph_data), rownames(edna_data))
cat("Number of morphological families:", nrow(morph_data), "\n")
cat("Number of eDNA families:", nrow(edna_data), "\n")
cat("Total number of families (union):", length(all_taxa), "\n")

# Fill missing families with zeros
morph_data_all <- morph_data %>%
  tibble::rownames_to_column("Family") %>%
  right_join(data.frame(Family = all_taxa), by = "Family") %>%
  replace(is.na(.), 0) %>%
  tibble::column_to_rownames("Family")

edna_data_all <- edna_data %>%
  tibble::rownames_to_column("Family") %>%
  right_join(data.frame(Family = all_taxa), by = "Family") %>%
  replace(is.na(.), 0) %>%
  tibble::column_to_rownames("Family")

##############################################
# Transpose function (sites as rows)
##############################################
preprocess_transposed_data <- function(df) {
  df_t <- as.data.frame(t(df))
  df_t$SiteID <- rownames(df_t)
  df_t <- df_t %>% select(SiteID, everything())
  return(df_t)
}

morph_transposed <- preprocess_transposed_data(morph_data_all)
edna_transposed  <- preprocess_transposed_data(edna_data_all)

cat("\nDimensions after transposition:\n")
cat("Morphological data:", paste(dim(morph_transposed), collapse = " x "), "\n")
cat("eDNA data:", paste(dim(edna_transposed), collapse = " x "), "\n")

##############################################
# Calculate Bray-Curtis similarity
##############################################
calculate_bc_similarity_transposed <- function(morph_df, edna_df, taxa_names) {
  results <- data.frame(
    SiteID = morph_df$SiteID,
    BrayCurtis_Similarity = NA_real_,
    Morph_Taxa_Count = NA_integer_,
    eDNA_Taxa_Count = NA_integer_,
    Shared_Taxa_Count = NA_integer_
  )
  
  for (i in seq_len(nrow(morph_df))) {
    site <- morph_df$SiteID[i]
    morph_vec <- as.numeric(morph_df[i, taxa_names])
    edna_vec  <- as.numeric(edna_df[i, taxa_names])
    
    morph_count <- sum(morph_vec > 0)
    edna_count  <- sum(edna_vec > 0)
    shared_count <- sum(morph_vec > 0 & edna_vec > 0)
    
    if (sum(morph_vec) == 0 & sum(edna_vec) == 0) {
      bc_similarity <- 1
    } else {
      community_matrix <- rbind(morph_vec, edna_vec)
      bc_similarity <- 1 - as.numeric(vegdist(community_matrix, method = "bray"))
    }
    if (is.na(bc_similarity)) bc_similarity <- 0
    
    results[i, ] <- list(site, bc_similarity, morph_count, edna_count, shared_count)
  }
  
  return(results)
}

cat("\n=== Calculating Bray-Curtis similarity (using all families) ===\n")

all_taxa_names <- colnames(morph_transposed)[-1]
bc_results <- calculate_bc_similarity_transposed(morph_transposed, edna_transposed, all_taxa_names)

##############################################
# Results summary and visualization
##############################################
cat("\n=== Calculation complete, generating summary statistics ===\n")

summary_stats <- bc_results %>%
  summarise(
    N = n(),
    Mean = mean(BrayCurtis_Similarity),
    SD = sd(BrayCurtis_Similarity),
    Median = median(BrayCurtis_Similarity),
    Min = min(BrayCurtis_Similarity),
    Max = max(BrayCurtis_Similarity),
    Q25 = quantile(BrayCurtis_Similarity, 0.25),
    Q75 = quantile(BrayCurtis_Similarity, 0.75)
  )
print(summary_stats)

# Plotting
p1 <- ggplot(bc_results, aes(x = BrayCurtis_Similarity)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = median(BrayCurtis_Similarity)), color = "red", linetype = "dashed") +
  labs(title = "Distribution of Bray-Curtis similarity between morphological and eDNA family-level composition (using all families)",
       x = "Bray-Curtis similarity", y = "Number of sites") +
  theme_minimal()

p2 <- ggplot(bc_results, aes(x = Shared_Taxa_Count, y = BrayCurtis_Similarity)) +
  geom_point(aes(size = Morph_Taxa_Count), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Relationship between similarity and number of shared families",
       x = "Number of shared families", y = "Bray-Curtis similarity", size = "Number of morphological families") +
  theme_minimal()

print(p1)
print(p2)

# Save results
##############################################
write.csv(bc_results, "BC_similarity_scores_allTaxa.csv", row.names = FALSE)

# Identify high-consistency sites (similarity > 0.5)
high_consistency_sites <- bc_results$SiteID[bc_results$BrayCurtis_Similarity > 0.5]
write.csv(data.frame(SiteID = high_consistency_sites),
          "high_consistency_sites_allTaxa.csv", row.names = FALSE)

cat("\n Analysis complete! Results saved:\n")
cat(" - BC_similarity_scores_allTaxa.csv (all sites)\n")
cat(" - high_consistency_sites_allTaxa.csv (high-consistency sites)\n")
##############################################
