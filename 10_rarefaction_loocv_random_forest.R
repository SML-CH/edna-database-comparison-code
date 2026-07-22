############################################################
## LOOCV with Rarefaction normalization (no data leakage)
## Unannotated OTUs predicting BMWP
############################################################

# Load required packages
library(e1071)
library(dplyr)
library(tidyr)
library(viridis)
library(ranger)
library(caret)
library(readxl)
library(ggplot2)
library(ggrepel)

##### Set working directory #####
setwd("E:/PhD/Database_comparison/eDNA_submission/Revision_results")

##### Read data #####
input_file <- "Unannotated-eDNA-BMWP.csv"

data <- read.csv(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

##### Extract BMWP #####
BMWP <- as.numeric(data[nrow(data), ])

##### Remove last row (BMWP), keep OTU count table #####
otu_counts <- data[-nrow(data), ]

##### Convert to matrix and ensure non-negative integer type #####
otu_counts <- as.matrix(otu_counts)
mode(otu_counts) <- "numeric"
otu_counts[is.na(otu_counts)] <- 0
otu_counts[otu_counts < 0] <- 0
otu_counts <- round(otu_counts)

##### Check data format: rows = OTUs, columns = samples #####
if (ncol(otu_counts) != length(BMWP)) {
  stop("Sample count mismatch, please check data format")
}

sample_ids <- colnames(otu_counts)

cat("Number of samples:", ncol(otu_counts), "\n")
cat("Number of OTUs:", nrow(otu_counts), "\n")

############################################################
## Function: Rarefaction for one sample
############################################################

rarefy_one_sample <- function(count_vector, depth) {
  total_reads <- sum(count_vector)
  
  if (total_reads < depth) {
    stop("Sample total reads below rarefaction depth: ", 
         total_reads, " < ", depth)
  }
  
  otu_index <- rep(seq_along(count_vector), count_vector)
  sampled_index <- sample(otu_index, size = depth, replace = FALSE)
  as.numeric(tabulate(sampled_index, nbins = length(count_vector)))
}

############################################################
## Function: Training-set Rarefaction normalization (no information leakage)
## Input:
## train_x: rows = OTUs, columns = training samples
## test_x : rows = OTUs, columns = test samples
##
## Output:
## train_rare: rows = samples, columns = OTUs, for ranger modeling
## test_rare : rows = samples, columns = OTUs, for ranger prediction
############################################################

rarefaction_normalize_train_test <- function(train_x, test_x, seed_offset = 0) {
  
  ##### 1. Remove OTUs with all zeros in training set #####
  keep_otus <- rowSums(train_x, na.rm = TRUE) > 0
  
  train_x <- train_x[keep_otus, , drop = FALSE]
  test_x  <- test_x[keep_otus, , drop = FALSE]
  
  ##### 2. Calculate rarefaction depth using training set (minimum sequencing depth) #####
  train_depths <- colSums(train_x)
  
  if (any(train_depths <= 0)) {
    stop("Training set contains samples with total reads = 0")
  }
  
  rarefaction_depth <- min(train_depths)
  
  ##### 3. Rarefy training set #####
  set.seed(1234 + seed_offset)
  
  train_rare_matrix <- apply(
    train_x,
    2,
    rarefy_one_sample,
    depth = rarefaction_depth
  )
  
  rownames(train_rare_matrix) <- rownames(train_x)
  colnames(train_rare_matrix) <- colnames(train_x)
  
  ##### 4. Rarefy test set with the same depth #####
  test_depths <- colSums(test_x)
  
  if (any(test_depths < rarefaction_depth)) {
    warning("Test sample depth (", min(test_depths), 
            ") less than training rarefaction depth (", rarefaction_depth, 
            "), using test set minimum depth")
    
    actual_depth <- min(test_depths)
    
    if (actual_depth <= 0) {
      stop("Test set has samples with zero depth, cannot rarefy")
    }
    
    set.seed(1234 + seed_offset + 1000)
    
    test_rare_matrix <- apply(
      test_x,
      2,
      rarefy_one_sample,
      depth = actual_depth
    )
    
    attr(test_rare_matrix, "rarefaction_depth") <- actual_depth
  } else {
    set.seed(1234 + seed_offset + 1000)
    
    test_rare_matrix <- apply(
      test_x,
      2,
      rarefy_one_sample,
      depth = rarefaction_depth
    )
    
    attr(test_rare_matrix, "rarefaction_depth") <- rarefaction_depth
  }
  
  rownames(test_rare_matrix) <- rownames(test_x)
  colnames(test_rare_matrix) <- colnames(test_x)
  
  ##### Transpose to rows = samples, columns = OTUs #####
  train_rare <- as.data.frame(t(train_rare_matrix), check.names = FALSE)
  test_rare <- as.data.frame(t(test_rare_matrix), check.names = FALSE)
  
  return(list(
    train_rare = train_rare,
    test_rare = test_rare,
    rarefaction_depth = rarefaction_depth
  ))
}

############################################################
## LOOCV main loop 
############################################################

set.seed(1234)

# Number of OTUs to retain for testing
top_n_values <- c(10, 20, 50, 100, 200, 500)
num_trees <- 500

# Save predictions for each top_n across all folds
all_predictions <- list()

# Save OTU importance for each fold
all_selected_otus <- list()

# Save rarefaction depth info for each fold
rarefaction_depth_info <- list()

n_sample <- ncol(otu_counts)

for (i in seq_len(n_sample)) {
  
  cat("LOOCV fold:", i, "/", n_sample, "\n")
  
  ##### 1. Split training and test sets (split on original data) #####
  train_index <- setdiff(seq_len(n_sample), i)
  test_index <- i
  
  train_raw <- otu_counts[, train_index, drop = FALSE]
  test_raw <- otu_counts[, test_index, drop = FALSE]
  
  train_y <- BMWP[train_index]
  test_y <- BMWP[test_index]
  
  ##### 2. Training-set Rarefaction normalization (no leakage!) #####
  rare_res <- rarefaction_normalize_train_test(train_raw, test_raw, seed_offset = i)
  
  train_x <- rare_res$train_rare  # rows = samples, columns = OTUs
  test_x <- rare_res$test_rare    # rows = samples, columns = OTUs
  
  # Record rarefaction depth
  rarefaction_depth_info[[i]] <- data.frame(
    Fold = i,
    Left_out_sample = sample_ids[i],
    Rarefaction_depth = rare_res$rarefaction_depth
  )
  
  ##### 3. Remove OTUs with zero variance in training set #####
  otu_var <- apply(train_x, 2, var, na.rm = TRUE)
  keep_var <- otu_var > 0 & !is.na(otu_var)
  
  train_x <- train_x[, keep_var, drop = FALSE]
  test_x <- test_x[, colnames(train_x), drop = FALSE]
  
  ##### If remaining OTUs are fewer than the minimum top_n, skip #####
  if (ncol(train_x) < min(top_n_values)) {
    cat("Warning: Fold", i, "has insufficient OTUs (", ncol(train_x), 
        ") for the minimum top_n. Skipping this fold.\n")
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
    seed = 1234 + i
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
    test_top <- test_x[, important_otus, drop = FALSE]
    
    mtry_top <- max(1, floor(ncol(train_top) / 3))
    
    rf_model <- ranger(
      x = train_top,
      y = train_y,
      num.trees = num_trees,
      importance = "permutation",
      mtry = mtry_top,
      min.node.size = 1,
      splitrule = "variance",
      seed = 4321 + i + top_n
    )
    
    pred <- predict(rf_model, data = test_top)$predictions
    
    all_predictions[[length(all_predictions) + 1]] <- data.frame(
      Normalization = "Rarefaction",
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
  "Unannotated-BMWP_Rarefaction_no_leakage_predictions.csv",
  row.names = FALSE
)

############################################################
## Save rarefaction depth info for each fold
############################################################

rarefaction_depth_df <- bind_rows(rarefaction_depth_info)

write.csv(
  rarefaction_depth_df,
  "Unannotated-BMWP_Rarefaction_no_leakage_rarefaction_depth.csv",
  row.names = FALSE
)

cat("\n=== Rarefaction depth summary ===\n")
summary(rarefaction_depth_df$Rarefaction_depth)

############################################################
## Save OTU importance for each fold
############################################################

importance_df <- bind_rows(all_selected_otus)

write.csv(
  importance_df,
  "Unannotated-BMWP_Rarefaction_no_leakage_LOOCV_OTU_importance.csv",
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
  ) %>%
  mutate(Normalization = "Rarefaction (no leakage)", .before = Top_n)

print(performance_df)

write.csv(
  performance_df,
  "Unannotated-BMWP_Rarefaction_no_leakage_feature_number_performance.csv",
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

p_feature_curve <- ggplot(performance_long, aes(x = Top_n, y = Value)) +
  geom_line(linewidth = 0.9, color = "darkgreen") +
  geom_point(size = 2.8, color = "darkgreen") +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = top_n_values) +
  labs(
    title = "LOOCV with within-training-set Rarefaction (no data leakage)",
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
  "Unannotated-BMWP_Rarefaction_no_leakage_feature_number_performance_curve.png",
  p_feature_curve,
  width = 7,
  height = 8,
  dpi = 300
)

############################################################
## Plot R2 curve separately
############################################################

p_r2_curve <- ggplot(performance_df, aes(x = Top_n, y = R2)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_point(size = 3, color = "darkgreen") +
  scale_x_continuous(breaks = top_n_values) +
  labs(
    title = "LOOCV with within-training-set Rarefaction (no data leakage)",
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
  "Unannotated-BMWP_Rarefaction_no_leakage_R2_feature_number_curve.png",
  p_r2_curve,
  width = 6,
  height = 4.5,
  dpi = 300
)
