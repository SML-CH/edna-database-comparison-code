script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}
source(file.path(script_dir, "_shared.R"), encoding = "UTF-8")

ensure_packages(c("vegan", "dplyr", "ggplot2"))

project_dir <- project_dir_from_script(script_dir)
data_dir <- data_dir_from_env(project_dir)
output_dir <- output_dir_from_env(project_dir, "bray_curtis_similarity")

morphology_file <- Sys.getenv(
  "MORPHOLOGY_FAMILY_FILE",
  unset = file.path(data_dir, "zhejiang_morphology_family_data.csv")
)

edna_file <- Sys.getenv(
  "EDNA_FAMILY_FILE",
  unset = file.path(data_dir, "NCBI_family_data.csv")
)

calculate_bray_curtis_similarity <- function(morphology_data, edna_data) {
  all_taxa <- union(rownames(morphology_data), rownames(edna_data))
  all_sites <- union(colnames(morphology_data), colnames(edna_data))

  morphology_aligned <- align_matrix(morphology_data, all_taxa, all_sites)
  edna_aligned <- align_matrix(edna_data, all_taxa, all_sites)

  results <- data.frame(
    SiteID = all_sites,
    BrayCurtis_Similarity = NA_real_,
    Morphology_Taxa_Count = NA_integer_,
    eDNA_Taxa_Count = NA_integer_,
    Shared_Taxa_Count = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(all_sites)) {
    site <- all_sites[[i]]
    morphology_vector <- as.numeric(morphology_aligned[, site])
    edna_vector <- as.numeric(edna_aligned[, site])

    morphology_count <- sum(morphology_vector > 0)
    edna_count <- sum(edna_vector > 0)
    shared_count <- sum(morphology_vector > 0 & edna_vector > 0)

    if (sum(morphology_vector) == 0 && sum(edna_vector) == 0) {
      similarity <- 1
    } else {
      community_matrix <- rbind(morphology_vector, edna_vector)
      similarity <- 1 - as.numeric(vegan::vegdist(community_matrix, method = "bray"))
    }

    if (is.na(similarity)) {
      similarity <- 0
    }

    results[i, ] <- list(site, similarity, morphology_count, edna_count, shared_count)
  }

  results
}

morphology_data <- read_numeric_matrix_csv(morphology_file)
edna_data <- read_numeric_matrix_csv(edna_file)

bc_results <- calculate_bray_curtis_similarity(morphology_data, edna_data)

summary_stats <- bc_results |>
  dplyr::summarise(
    N = dplyr::n(),
    Mean = mean(BrayCurtis_Similarity),
    SD = sd(BrayCurtis_Similarity),
    Median = median(BrayCurtis_Similarity),
    Min = min(BrayCurtis_Similarity),
    Max = max(BrayCurtis_Similarity),
    Q25 = quantile(BrayCurtis_Similarity, 0.25),
    Q75 = quantile(BrayCurtis_Similarity, 0.75)
  )

histogram_plot <- ggplot2::ggplot(bc_results, ggplot2::aes(x = BrayCurtis_Similarity)) +
  ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
  ggplot2::geom_vline(
    ggplot2::aes(xintercept = median(BrayCurtis_Similarity)),
    color = "red",
    linetype = "dashed"
  ) +
  ggplot2::labs(
    title = "Family-level Bray-Curtis similarity between morphology and eDNA",
    x = "Bray-Curtis similarity",
    y = "Number of sites"
  ) +
  ggplot2::theme_minimal()

shared_taxa_plot <- ggplot2::ggplot(
  bc_results,
  ggplot2::aes(x = Shared_Taxa_Count, y = BrayCurtis_Similarity)
) +
  ggplot2::geom_point(ggplot2::aes(size = Morphology_Taxa_Count), alpha = 0.6) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red") +
  ggplot2::labs(
    title = "Similarity and shared family count",
    x = "Number of shared families",
    y = "Bray-Curtis similarity",
    size = "Morphology families"
  ) +
  ggplot2::theme_minimal()

write.csv(
  bc_results,
  file.path(output_dir, "bray_curtis_similarity_scores_all_taxa.csv"),
  row.names = FALSE
)

write.csv(
  summary_stats,
  file.path(output_dir, "bray_curtis_similarity_summary.csv"),
  row.names = FALSE
)

ggplot2::ggsave(
  file.path(output_dir, "bray_curtis_similarity_histogram.pdf"),
  histogram_plot,
  width = 7,
  height = 5
)

ggplot2::ggsave(
  file.path(output_dir, "bray_curtis_similarity_shared_taxa.pdf"),
  shared_taxa_plot,
  width = 7,
  height = 5
)

threshold_text <- Sys.getenv("BC_HIGH_CONSISTENCY_THRESHOLD", unset = "")
if (nzchar(threshold_text)) {
  high_consistency_threshold <- as.numeric(threshold_text)

  if (is.na(high_consistency_threshold)) {
    stop("BC_HIGH_CONSISTENCY_THRESHOLD must be numeric.", call. = FALSE)
  }

  high_consistency_sites <- bc_results |>
    dplyr::filter(BrayCurtis_Similarity >= high_consistency_threshold)

  write.csv(
    high_consistency_sites,
    file.path(output_dir, "high_consistency_sites_all_taxa.csv"),
    row.names = FALSE
  )
}

print(summary_stats)
