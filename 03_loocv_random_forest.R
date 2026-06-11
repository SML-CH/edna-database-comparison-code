script_dir <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}
source(file.path(script_dir, "_shared.R"), encoding = "UTF-8")

ensure_packages(c("dplyr", "ranger", "metagenomeSeq", "Biobase"))

project_dir <- project_dir_from_script(script_dir)
data_dir <- data_dir_from_env(project_dir)
output_dir <- output_dir_from_env(project_dir, "loocv_random_forest_wqi")

input_file <- Sys.getenv("WQI_OTU_FILE", unset = file.path(data_dir, "CMBL-WQI.csv"))
response_row <- Sys.getenv("WQI_RESPONSE_ROW", unset = "WQI")
top_n <- as.integer(Sys.getenv("TOP_N_OTUS", unset = "200"))
num_trees <- as.integer(Sys.getenv("RF_NUM_TREES", unset = "500"))
random_seed <- as.integer(Sys.getenv("RF_RANDOM_SEED", unset = "1234"))

y <- as.numeric(data[response_row, ])
otu_raw <- data[rownames(data) != response_row, , drop = FALSE]
otu_raw <- as.data.frame(otu_raw, check.names = FALSE)
otu_raw[] <- lapply(otu_raw, function(x) as.numeric(as.character(x)))

sample_ids <- colnames(otu_raw)

if (length(y) != ncol(otu_raw)) {
  stop("The response length does not match the number of samples.", call. = FALSE)
}

if (any(is.na(y))) {
  stop("The response row contains NA values.", call. = FALSE)
}

if (any(is.na(otu_raw))) {
  stop("The OTU matrix contains NA values.", call. = FALSE)
}

css_normalize_train_test <- function(train_x, test_x) {
  keep_otus <- rowSums(train_x, na.rm = TRUE) > 0

  train_x <- train_x[keep_otus, , drop = FALSE]
  test_x <- test_x[keep_otus, , drop = FALSE]

  train_counts <- as.matrix(train_x)
  train_pheno <- Biobase::AnnotatedDataFrame(
    data.frame(row.names = colnames(train_counts))
  )

  mr_train <- metagenomeSeq::newMRexperiment(
    counts = train_counts,
    phenoData = train_pheno
  )

  p_train <- metagenomeSeq::cumNormStatFast(mr_train)
  mr_train <- metagenomeSeq::cumNorm(mr_train, p = p_train)
  train_css_matrix <- metagenomeSeq::MRcounts(mr_train, norm = TRUE, log = FALSE)

  test_counts <- as.matrix(test_x)
  test_pheno <- Biobase::AnnotatedDataFrame(
    data.frame(row.names = colnames(test_counts))
  )

  mr_test <- metagenomeSeq::newMRexperiment(
    counts = test_counts,
    phenoData = test_pheno
  )

  mr_test <- metagenomeSeq::cumNorm(mr_test, p = p_train)
  test_css_matrix <- metagenomeSeq::MRcounts(mr_test, norm = TRUE, log = FALSE)

  list(
    train_css = as.data.frame(t(train_css_matrix), check.names = FALSE),
    test_css = as.data.frame(t(test_css_matrix), check.names = FALSE)
  )
}

set.seed(random_seed)

predictions_by_fold <- list()
importance_by_fold <- list()
n_sample <- ncol(otu_raw)

for (i in seq_len(n_sample)) {
  message("LOOCV fold: ", i, "/", n_sample)

  train_index <- setdiff(seq_len(n_sample), i)
  test_index <- i

  train_x_raw <- otu_raw[, train_index, drop = FALSE]
  test_x_raw <- otu_raw[, test_index, drop = FALSE]

  train_y <- y[train_index]
  test_y <- y[test_index]

  css_result <- css_normalize_train_test(train_x_raw, test_x_raw)

  train_x <- css_result$train_css
  test_x <- css_result$test_css

  otu_variance <- apply(train_x, 2, var, na.rm = TRUE)
  keep_variance <- otu_variance > 0 & !is.na(otu_variance)

  train_x <- train_x[, keep_variance, drop = FALSE]
  test_x <- test_x[, keep_variance, drop = FALSE]

  if (ncol(train_x) < 1) {
    stop("No OTUs are available for model training in fold ", i, ".", call. = FALSE)
  }

  mtry_all <- max(1, floor(ncol(train_x) / 3))

  rf_importance <- ranger::ranger(
    x = train_x,
    y = train_y,
    num.trees = num_trees,
    importance = "permutation",
    mtry = mtry_all,
    min.node.size = 1,
    splitrule = "variance",
    seed = random_seed
  )

  importance_scores <- ranger::importance(rf_importance)

  importance_df <- data.frame(
    Fold = i,
    Left_out_sample = sample_ids[i],
    OTU = names(importance_scores),
    Importance = as.numeric(importance_scores),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(Importance))

  n_select <- min(top_n, nrow(importance_df))
  selected_otus <- importance_df$OTU[seq_len(n_select)]
  importance_df$Selected_in_this_fold <- importance_df$OTU %in% selected_otus

  importance_by_fold[[i]] <- importance_df

  train_x_top <- train_x[, selected_otus, drop = FALSE]
  test_x_top <- test_x[, selected_otus, drop = FALSE]

  mtry_top <- max(1, floor(ncol(train_x_top) / 3))

  rf_top <- ranger::ranger(
    x = train_x_top,
    y = train_y,
    num.trees = num_trees,
    importance = "permutation",
    mtry = mtry_top,
    min.node.size = 1,
    splitrule = "variance",
    seed = random_seed
  )

  pred_value <- predict(rf_top, data = test_x_top)$predictions

  predictions_by_fold[[i]] <- data.frame(
    SampleID = sample_ids[i],
    Observed_WQI = test_y,
    Predicted_WQI = pred_value,
    Fold = i,
    Number_selected_OTUs = n_select
  )
}

prediction_results <- dplyr::bind_rows(predictions_by_fold)
r2_value <- cor(prediction_results$Observed_WQI, prediction_results$Predicted_WQI)^2

performance <- data.frame(
  Model = "Random forest",
  Cross_validation = "LOOCV",
  Normalization = "CSS within training set",
  Importance = "Permutation importance",
  Feature_selection = paste0("Top ", top_n, " OTUs selected within each training set"),
  R2 = r2_value
)

importance_all_folds <- dplyr::bind_rows(importance_by_fold)

otu_selection_summary <- importance_all_folds |>
  dplyr::group_by(OTU) |>
  dplyr::summarise(
    Selected_frequency = sum(Selected_in_this_fold),
    Selected_ratio = mean(Selected_in_this_fold),
    Mean_importance = mean(Importance, na.rm = TRUE),
    Median_importance = median(Importance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(Selected_frequency), dplyr::desc(Mean_importance))

write.csv(
  prediction_results,
  file.path(output_dir, "cmbl_wqi_loocv_predictions_css_top_otus_permutation.csv"),
  row.names = FALSE
)

write.csv(
  performance,
  file.path(output_dir, "cmbl_wqi_loocv_model_performance_css_top_otus_permutation.csv"),
  row.names = FALSE
)

write.csv(
  importance_all_folds,
  file.path(output_dir, "cmbl_wqi_loocv_fold_permutation_importance.csv"),
  row.names = FALSE
)

write.csv(
  otu_selection_summary,
  file.path(output_dir, "cmbl_wqi_otu_selection_frequency_permutation.csv"),
  row.names = FALSE
)

print(performance)

