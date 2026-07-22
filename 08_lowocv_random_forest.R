############################################################
## Leave-one-watershed-out validation
## raw OTU counts -> leave one Watershed out
## -> training-set CSS normalization
## -> training-set OTU importance ranking
## -> automatic feature-number selection by training OOB RMSE
## -> prediction on the left-out watershed
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(metagenomeSeq)
library(Biobase)

setwd("E:/PhD/Database_comparison/China_freshwater_benthic_database")

##### load data#####
raw_data <- read.csv(
  "Taxonomy-free_WQI_Watershed.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta_rows <- c("WQI", "Watershed")
otu_rows <- setdiff(rownames(raw_data), meta_rows)

otu_data <- as.data.frame(t(raw_data[otu_rows, , drop = FALSE]))
otu_data <- as.data.frame(lapply(otu_data, as.numeric))
rownames(otu_data) <- colnames(raw_data)

meta_data <- as.data.frame(t(raw_data[meta_rows, , drop = FALSE]))
rownames(meta_data) <- colnames(raw_data)

data <- cbind(otu_data, meta_data)

data$WQI <- as.numeric(data$WQI)
data$Watershed <- as.factor(data$Watershed)
##### set  parameter#####
set.seed(1234)

dataset_tag <- "Taxonomy-free-WQI"

y_col <- "WQI"
group_col <- "Watershed"

feature_numbers <- c(20, 50, 100, 200, 500)

meta_cols <- intersect(c(y_col, group_col, "Year"), colnames(data))
otu_cols <- setdiff(colnames(data), meta_cols)

##### extract raw OTU counts #####
otu_raw_all <- data[, otu_cols, drop = FALSE]
otu_raw_all <- as.data.frame(lapply(otu_raw_all, as.numeric))
rownames(otu_raw_all) <- rownames(data)

##### CSS #####
css_normalize_train_test <- function(train_counts, test_counts) {
  
  # train_counts/test_counts: samples x OTUs
  
  train_counts <- as.matrix(train_counts)
  test_counts  <- as.matrix(test_counts)
  
  storage.mode(train_counts) <- "numeric"
  storage.mode(test_counts)  <- "numeric"
  
  # metagenomeSeq 要求 OTUs x samples
  train_mat <- t(train_counts)
  test_mat  <- t(test_counts)
  
  pheno_train <- AnnotatedDataFrame(
    data.frame(row.names = colnames(train_mat))
  )
  
  mr_train <- newMRexperiment(
    counts = train_mat,
    phenoData = pheno_train
  )
  
  # calculate CSS percentile
  p_train <- cumNormStatFast(mr_train)
  
  mr_train <- cumNorm(mr_train, p = p_train)
  
  train_css <- MRcounts(
    mr_train,
    norm = TRUE,
    log = FALSE
  )
  
  train_css <- as.data.frame(t(train_css))
  
 
  pheno_test <- AnnotatedDataFrame(
    data.frame(row.names = colnames(test_mat))
  )
  
  mr_test <- newMRexperiment(
    counts = test_mat,
    phenoData = pheno_test
  )
  
  mr_test <- cumNorm(mr_test, p = p_train)
  
  test_css <- MRcounts(
    mr_test,
    norm = TRUE,
    log = FALSE
  )
  
  test_css <- as.data.frame(t(test_css))
  
  test_css <- test_css[, colnames(train_css), drop = FALSE]
  
  return(
    list(
      train_css = train_css,
      test_css = test_css,
      p_train = p_train
    )
  )
}

############################################################
## Leave-one-watershed-out validation
## raw OTU counts -> leave one Watershed out
## -> training-set CSS normalization
## -> training-set OTU importance ranking
## -> evaluate each feature number on left-out watershed
## -> only save prediction performance for each feature number
############################################################

watersheds <- levels(data[[group_col]])

all_predictions <- list()

for (ws in watersheds) {
  
  cat("\n==============================\n")
  cat("leave watershed:", ws, "\n")
  cat("==============================\n")
  
  train_id <- data[[group_col]] != ws
  test_id  <- data[[group_col]] == ws
  
  train_meta <- data[train_id, , drop = FALSE]
  test_meta  <- data[test_id, , drop = FALSE]
  
  train_raw <- otu_raw_all[train_id, , drop = FALSE]
  test_raw  <- otu_raw_all[test_id, , drop = FALSE]
  
  ##### 1. filter OTU #####
  keep_nonzero <- colSums(train_raw, na.rm = TRUE) > 0
  train_raw <- train_raw[, keep_nonzero, drop = FALSE]
  test_raw  <- test_raw[, colnames(train_raw), drop = FALSE]
  
  keep_var <- apply(train_raw, 2, var, na.rm = TRUE) > 0
  train_raw <- train_raw[, keep_var, drop = FALSE]
  test_raw  <- test_raw[, colnames(train_raw), drop = FALSE]
  
  cat("train sites:", nrow(train_raw), "\n")
  cat("test sites:", nrow(test_raw), "\n")
  cat(" OTU number:", ncol(train_raw), "\n")
  
  ##### 2. train CSS #####
  css_out <- css_normalize_train_test(
    train_counts = train_raw,
    test_counts  = test_raw
  )
  
  x_train <- css_out$train_css
  x_test  <- css_out$test_css
  y_train <- train_meta[[y_col]]
  
  cat("train CSS percentile p:", css_out$p_train, "\n")
  
  ##### 3. train OTU importance  #####
  rf_importance <- ranger(
    x = x_train,
    y = y_train,
    num.trees = 500,
    importance = "permutation",
    mtry = max(1, floor(ncol(x_train) / 3)),
    min.node.size = 1,
    seed = 1234
  )
  
  importance_scores <- importance(rf_importance)
  
  importance_df <- data.frame(
    OTU = names(importance_scores),
    Importance = as.numeric(importance_scores),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Importance))
  
  ##### 4. Model and predict held-out watersheds separately for each feature count #####
  candidate_numbers <- unique(c(
    feature_numbers[feature_numbers < ncol(x_train)],
    ncol(x_train)
  ))
  
  for (top_n in candidate_numbers) {
    
    feature_label <- ifelse(
      top_n == ncol(x_train),
      "All",
      as.character(top_n)
    )
    
    cat("current feature count:", feature_label, "\n")
    
    selected_otus <- importance_df$OTU[1:top_n]
    
    x_train_top <- x_train[, selected_otus, drop = FALSE]
    x_test_top  <- x_test[, selected_otus, drop = FALSE]
    
    rf_model <- ranger(
      x = x_train_top,
      y = y_train,
      num.trees = 500,
      importance = "none",
      mtry = max(1, floor(ncol(x_train_top) / 3)),
      min.node.size = 1,
      seed = 4321 + top_n
    )
    
    pred <- predict(rf_model, data = x_test_top)$predictions
    
    all_predictions[[length(all_predictions) + 1]] <- data.frame(
      Feature_number = feature_label,
      Actual_feature_number = top_n,
      Left_out_watershed = ws,
      SampleID = rownames(test_meta),
      Observed_WQI = test_meta[[y_col]],
      Predicted_WQI = pred,
      OOB_RMSE_training = sqrt(rf_model$prediction.error),
      CSS_p_train = css_out$p_train,
      stringsAsFactors = FALSE
    )
  }
}

############################################################
############################################################

prediction_df <- bind_rows(all_predictions)

prediction_df$Feature_number <- factor(
  prediction_df$Feature_number,
  levels = c("20", "50", "100", "200", "500", "All")
)

############################################################
############################################################

performance_each_feature <- prediction_df %>%
  group_by(Feature_number, Actual_feature_number) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean((Observed_WQI - Predicted_WQI)^2, na.rm = TRUE)),
    MAE = mean(abs(Observed_WQI - Predicted_WQI), na.rm = TRUE),
    R2 = cor(Observed_WQI, Predicted_WQI, use = "complete.obs")^2,
    Mean_training_OOB_RMSE = mean(OOB_RMSE_training, na.rm = TRUE),
    SD_training_OOB_RMSE = sd(OOB_RMSE_training, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Actual_feature_number)

print(performance_each_feature)

write.csv(
  performance_each_feature,
  paste0(dataset_tag, "_LOWO_each-feature-number_performance.csv"),
  row.names = FALSE
)


