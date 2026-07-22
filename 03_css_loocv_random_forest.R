############################################################
## CMBL BMWP prediction with LOOCV and CSS normalization
## Sensitivity analysis of feature number (no data leakage)
############################################################

# Load required packages
library(e1071)
library(dplyr)
library(viridis)
library(ranger)
library(caret)
library(readxl)
library(metagenomeSeq)
library(ggplot2)
library(ggrepel)
library(Biobase)
library(tidyr)

##### Set working directory #####
setwd("E:/PhD/Database_comparison/China_freshwater_benthic_database")

##### Read data #####
data <- read.csv(
  "CMBL_BMWP.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

##### Extract BMWP #####
BMWP <- as.numeric(data[nrow(data), ])

##### Remove last row (BMWP), keep OTU count table #####
otu_counts <- data[-nrow(data), ]

##### Convert to matrix and ensure numeric type #####
otu_counts <- as.matrix(otu_counts)
mode(otu_counts) <- "numeric"

##### Transpose to OTU rows x sample columns (consistent with code 1) #####
otu_raw <- otu_counts  # Now OTU rows x sample columns

##### Sample IDs #####
sample_ids <- colnames(otu_raw)

############################################################
## Function: training-set CSS normalization (no information leakage)
## Input:
## train_x: OTU rows x training samples
## test_x : OTU rows x test samples
##
## Output:
## train_css: sample rows x OTU columns, for ranger modeling
## test_css : sample rows x OTU columns, for ranger prediction
############################################################

css_normalize_train_test <- function(train_x, test_x) {
  
  ##### Remove OTUs with all zeros in training set #####
  keep_otus <- rowSums(train_x, na.rm = TRUE) > 0
  
  train_x <- train_x[keep_otus, , drop = FALSE]
  test_x  <- test_x[keep_otus, , drop = FALSE]
  
  ##### 1. Calculate CSS normalization parameter p using only training set #####
  train_counts <- as.matrix(train_x)
  
  train_pheno <- AnnotatedDataFrame(
    data.frame(row.names = colnames(train_counts))
  )
  
  mr_train <- newMRexperiment(
    counts = train_counts,
    phenoData = train_pheno
  )
  
  p_train <- cumNormStatFast(mr_train)
  
  ##### 2. Normalize training set using training-set p value #####
  mr_train <- cumNorm(mr_train, p = p_train)
  train_css_matrix <- MRcounts(mr_train, norm = TRUE, log = FALSE)
  
  ##### 3. Normalize test set using the same p value #####
  test_counts <- as.matrix(test_x)
  
  test_pheno <- AnnotatedDataFrame(
    data.frame(row.names = colnames(test_counts))
  )
  
  mr_test <- newMRexperiment(
    counts = test_counts,
    phenoData = test_pheno
  )
  
  mr_test <- cumNorm(mr_test, p = p_train)
  test_css_matrix <- MRcounts(mr_test, norm = TRUE, log = FALSE)
  
  ##### metagenomeSeq output: OTU rows x sample columns
  ##### ranger modeling requires: sample rows x OTU columns
  train_css <- as.data.frame(t(train_css_matrix), check.names = FALSE)
  test_css  <- as.data.frame(t(test_css_matrix), check.names = FALSE)
  
  return(list(
    train_css = train_css,
    test_css = test_css
  ))
}

############################################################
## LOOCV main loop 
############################################################

set.seed(1234)
total_otus <- nrow(otu_counts)
ALL <- total_otus  # ALL represents all features
# Number of OTUs to retain for testing
top_n_values <- c(20, 50, 100, 200, 500, ALL)
num_trees <- 500

# Save predictions for each top_n across all folds
all_predictions <- list()

# Save OTU importance for each fold
all_selected_otus <- list()

n_sample <- ncol(otu_raw)  # Number of samples (columns)

for (i in seq_len(n_sample)) {
  
  cat("LOOCV fold:", i, "/", n_sample, "\n")
  
  ##### 1. Leave out the i-th site as test set #####
  train_index <- setdiff(seq_len(n_sample), i)
  test_index  <- i
  
  ##### Note: data is OTU rows x sample columns, split by columns #####
  train_x_raw <- otu_raw[, train_index, drop = FALSE]
  test_x_raw  <- otu_raw[, test_index, drop = FALSE]
  
  train_y <- BMWP[train_index]
  test_y  <- BMWP[test_index]
  
  ##### 2. Training-set CSS normalization (no leakage!) #####
  css_res <- css_normalize_train_test(train_x_raw, test_x_raw)
  
  train_x <- css_res$train_css  # Sample rows x OTU columns
  test_x  <- css_res$test_css   # Sample rows x OTU columns
  
  ##### 3. Remove OTUs with zero variance in training set #####
  otu_var <- apply(train_x, 2, var, na.rm = TRUE)
  keep_var <- otu_var > 0 & !is.na(otu_var)
  
  train_x <- train_x[, keep_var, drop = FALSE]
  test_x  <- test_x[, keep_var, drop = FALSE]
  
  ##### If remaining OTUs are fewer than the minimum top_n, skip or warn #####
  if (ncol(train_x) < min(top_n_values)) {
    cat("Warning: Fold", i, "has insufficient OTUs (", ncol(train_x), ") for the minimum top_n. Skipping this fold.\n")
    next
  }
  
  ##########################################################
  ## Step 1: Calculate permutation importance using only training set
  ##########################################################
  
  mtry_all <- max(1, floor(ncol(train_x) / 3))
  
  rf_importance <- ranger(
    x = train_x,
    y = train_y,
    num.trees = num_trees,
    importance = "permutation",
    mtry = mtry_all,
    min.node.size = 1,
    splitrule = "variance",
    seed = 1234
  )
  
  importance_scores <- importance(rf_importance)
  
  importance_df_fold <- data.frame(
    Fold = i,
    TestSample = sample_ids[i],
    OTU = names(importance_scores),
    Importance = as.numeric(importance_scores)
  ) %>%
    arrange(desc(Importance))
  
  all_selected_otus[[i]] <- importance_df_fold
  
  ##########################################################
  ## Step 2: Model and predict for each top_n
  ##########################################################
  
  for (top_n in top_n_values) {
    
    top_k <- min(top_n, nrow(importance_df_fold))
    
    if (top_k == 0) {
      cat("Warning: Fold", i, ", top_n =", top_n, "has no selectable OTUs\n")
      next
    }
    
    important_otus <- importance_df_fold$OTU[1:top_k]
    
    train_top <- train_x[, important_otus, drop = FALSE]
    test_top  <- test_x[, important_otus, drop = FALSE]
    
    mtry_top <- max(1, floor(ncol(train_top) / 3))
    
    rf_model <- ranger(
      x = train_top,
      y = train_y,
      num.trees = num_trees,
      importance = "permutation",
      mtry = mtry_top,
      min.node.size = 1,
      splitrule = "variance",
      seed = 1234
    )
    
    pred <- predict(rf_model, data = test_top)$predictions
    
    all_predictions[[length(all_predictions) + 1]] <- data.frame(
      Top_n = top_n,
      Fold = i,
      SampleID = sample_ids[i],
      Observed_BMWP = test_y,
      Predicted_BMWP = pred
    )
  }
}

############################################################
## Compile all prediction results
############################################################

prediction_df <- bind_rows(all_predictions)

write.csv(
  prediction_df,
  "CMBL_BMWP_CSS_no_leakage_feature_number_predictions.csv",
  row.names = FALSE
)

############################################################
## Calculate model performance for different feature numbers
############################################################

performance_df <- prediction_df %>%
  group_by(Top_n) %>%
  summarise(
    RMSE = sqrt(mean((Observed_BMWP - Predicted_BMWP)^2)),
    MAE = mean(abs(Observed_BMWP - Predicted_BMWP)),
    R2 = cor(Observed_BMWP, Predicted_BMWP)^2,
    n_folds = n(),
    .groups = "drop"
  )

print(performance_df)

write.csv(
  performance_df,
  "CMBL_BMWP_CSS_no_leakage_feature_number_performance.csv",
  row.names = FALSE
)

############################################################
## Plot feature number vs model performance curves
############################################################

performance_long <- performance_df %>%
  pivot_longer(
    cols = c(RMSE, MAE, R2),
    names_to = "Metric",
    values_to = "Value"
  )

p_feature_curve <- ggplot(performance_long %>% filter(Top_n != ALL), aes(x = Top_n, y = Value)) +
  geom_line(linewidth = 0.9, color = "steelblue") +
  geom_point(size = 2.8, color = "steelblue") +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = top_n_values) +
  labs(
    title = "LOOCV with training-set CSS normalization (no data leakage)",
    x = "Number of retained OTUs",
    y = "Model performance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 11)
  )

print(p_feature_curve)

ggsave(
  "CMBL_BMWP_CSS_no_leakage_feature_number_performance_curve.png",
  p_feature_curve,
  width = 7,
  height = 8,
  dpi = 300
)

############################################################
## Plot R2 curve separately
############################################################

p_r2_curve <- ggplot(performance_df %>% filter(Top_n != ALL), aes(x = Top_n, y = R2)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  scale_x_continuous(breaks = top_n_values) +
  labs(
    title = "LOOCV with training-set CSS normalization (no data leakage)",
    x = "Number of retained OTUs",
    y = expression(R^2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 11)
  )

print(p_r2_curve)

ggsave(
  "CMBL_BMWP_CSS_no_leakage_R2_feature_number_curve.png",
  p_r2_curve,
  width = 6,
  height = 4.5,
  dpi = 300
)

############################################################
## Optional: Save importance for each fold
############################################################

importance_all_folds <- bind_rows(all_selected_otus)

write.csv(
  importance_all_folds,
  "CMBL_BMWP_CSS_no_leakage_fold_OTU_importance.csv",
  row.names = FALSE
)

cat("\n=== No-leakage version of feature number sensitivity analysis completed ===\n")
print(performance_df)

############################################################
## Optional: Compare with leakage version (if previously run)
############################################################

# if (file.exists("CMBL_BMWP_CSS_feature_number_performance.csv")) {
#   leaked_perf <- read.csv("CMBL_BMWP_CSS_feature_number_performance.csv")
#   leaked_perf$Version <- "With leakage"
#   no_leak_perf <- performance_df
#   no_leak_perf$Version <- "No leakage"
#   
#   compare_df <- rbind(
#     leaked_perf[, c("Top_n", "RMSE", "MAE", "R2", "Version")],
#     no_leak_perf[, c("Top_n", "RMSE", "MAE", "R2", "Version")]
#   )
#   
#   p_compare <- ggplot(compare_df, aes(x = Top_n, y = R2, color = Version)) +
#     geom_line(linewidth = 1) +
#     geom_point(size = 2.5) +
#     scale_x_continuous(breaks = top_n_values) +
#     labs(x = "Number of retained OTUs", y = expression(R^2)) +
#     theme_classic(base_size = 13) +
#     theme(axis.text = element_text(color = "black"))
#   
#   ggsave("CMBL_BMWP_CSS_leakage_comparison_R2_curve.png", p_compare, width = 6, height = 4.5, dpi = 300)
# }
