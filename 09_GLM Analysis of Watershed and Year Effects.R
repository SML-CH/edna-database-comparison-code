# ================================================================
# Watershed and year variance analysis (excluding Qiantang River; no eta squared; no result files)
# ================================================================
setwd("E:/PhD/Database_comparison/eDNA_revision_results")
# Install missing packages if needed:
required_packages <- c("readr", "dplyr", "purrr", "tidyr", "stringr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
})

# Set working directory to script location
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(script_path)) return(dirname(normalizePath(script_path)))
  }
  getwd()
}

setwd(get_script_dir())

# Six input files
input_files <- c(
  "Unannotated-WQI-Watershed.csv",
  "CMBL-WQI-Watershed.csv",
  "NCBI-WQI-Watershed.csv",
  "Unannotated-BMWP-Watershed.csv",
  "CMBL-BMWP-Watershed.csv",
  "NCBI-BMWP-Watershed.csv"
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Files not found:\n", paste(missing_files, collapse = "\n"))
}

# Extract R-squared
get_r2 <- function(model) {
  unname(summary(model)$r.squared)
}

# Compare reduced and full models, return delta R2 and p value
compare_models <- function(reduced_model, full_model) {
  added_df <- df.residual(reduced_model) - df.residual(full_model)
  
  if (!is.finite(added_df) || added_df <= 0) {
    return(tibble(delta_R2 = NA_real_, p_value = NA_real_))
  }
  
  delta_r2 <- get_r2(full_model) - get_r2(reduced_model)
  if (delta_r2 > -1e-12) delta_r2 <- max(delta_r2, 0)
  
  model_test <- anova(reduced_model, full_model)
  
  tibble(
    delta_R2 = delta_r2,
    p_value = model_test[["Pr(>F)"]][2]
  )
}

analyse_one_file <- function(file_name) {
  raw_data <- read_csv(
    file_name,
    show_col_types = FALSE,
    locale = locale(encoding = "UTF-8")
  )
  
  index_name <- if (str_detect(file_name, "WQI")) "WQI" else "BMWP"
  observed_column <- paste0("Observed_", index_name)
  predicted_column <- paste0("Predicted_", index_name)
  
  required_columns <- c(
    "SampleID", observed_column, predicted_column, "watershed", "year"
  )
  absent_columns <- setdiff(required_columns, names(raw_data))
  if (length(absent_columns) > 0) {
    stop(file_name, " missing columns: ", paste(absent_columns, collapse = ", "))
  }
  
  analysis_data <- raw_data %>%
    transmute(
      SampleID = as.character(.data$SampleID),
      Observed = as.numeric(.data[[observed_column]]),
      Predicted = as.numeric(.data[[predicted_column]]),
      watershed = as.character(.data$watershed),
      year = as.character(.data$year)
    ) %>%
    drop_na(Observed, Predicted, watershed, year) %>%
    filter(!str_detect(
      watershed,
      regex("Qiantang|Qiantangjiang|钱塘", ignore_case = TRUE)
    )) %>%
    mutate(
      watershed = droplevels(factor(watershed)),
      year = droplevels(factor(year))
    )
  
  if (nrow(analysis_data) < 5) {
    stop(file_name, " insufficient valid samples after excluding Qiantang River.")
  }
  
  # Observed model: Observed ~ watershed + year
  observed_full <- lm(
    Observed ~ watershed + year,
    data = analysis_data,
    na.action = na.fail
  )
  observed_without_watershed <- lm(
    Observed ~ year,
    data = analysis_data,
    na.action = na.fail
  )
  observed_without_year <- lm(
    Observed ~ watershed,
    data = analysis_data,
    na.action = na.fail
  )
  
  observed_watershed <- compare_models(
    observed_without_watershed,
    observed_full
  )
  observed_year <- compare_models(
    observed_without_year,
    observed_full
  )
  
  observed_null <- lm(Observed ~ 1, data = analysis_data)
  observed_total_p <- anova(observed_null, observed_full)[["Pr(>F)"]][2]
  
  # Prediction model: Predicted ~ Observed + watershed + year
  prediction_full <- lm(
    Predicted ~ Observed + watershed + year,
    data = analysis_data,
    na.action = na.fail
  )
  prediction_without_watershed <- lm(
    Predicted ~ Observed + year,
    data = analysis_data,
    na.action = na.fail
  )
  prediction_without_year <- lm(
    Predicted ~ Observed + watershed,
    data = analysis_data,
    na.action = na.fail
  )
  
  prediction_watershed <- compare_models(
    prediction_without_watershed,
    prediction_full
  )
  prediction_year <- compare_models(
    prediction_without_year,
    prediction_full
  )
  
  prediction_null <- lm(Predicted ~ 1, data = analysis_data)
  prediction_total_p <- anova(
    prediction_null,
    prediction_full
  )[["Pr(>F)"]][2]
  
  dataset_name <- str_remove(basename(file_name), "-Watershed\\.csv$")
  
  tibble(
    Dataset = dataset_name,
    Index = index_name,
    n_after_excluding_Qiantang = nrow(analysis_data),
    `Observed: total R2` = get_r2(observed_full),
    `Observed: total model p` = observed_total_p,
    `Observed: watershed unique ΔR2` = observed_watershed$delta_R2,
    `Observed: watershed p` = observed_watershed$p_value,
    `Observed: year unique ΔR2` = observed_year$delta_R2,
    `Observed: year p` = observed_year$p_value,
    `Prediction: total R2` = get_r2(prediction_full),
    `Prediction: total model p` = prediction_total_p,
    `Prediction: watershed extra ΔR2 after Observed` =
      prediction_watershed$delta_R2,
    `Prediction: watershed p` = prediction_watershed$p_value,
    `Prediction: year extra ΔR2 after Observed` = prediction_year$delta_R2,
    `Prediction: year p` = prediction_year$p_value
  )
}

# Run analysis on six datasets, display results in console only
result <- map_dfr(input_files, analyse_one_file)

options(tibble.width = Inf, pillar.sigfig = 6)
print(result, n = Inf, width = Inf)

# Uncomment to view in RStudio data viewer:
# View(result)
