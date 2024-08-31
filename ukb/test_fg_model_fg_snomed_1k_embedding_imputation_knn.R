library(glmnet)
library(dplyr)
library(survival)
library(RANN)
library(data.table)
library(Matrix)
source("snomed_embedding_functions.R")
# If embedding_table is defined already, use the existing value
if (!exists("embedding_table")) {
  embedding_table <- load_embedding_table()
}

model_results_single_fold <- function(fold_number) {
  model_files_path <- "../models/1k_snomed/"

  lasso_coefs <- read.csv(paste0(model_files_path, "1k_snomed_cox_lasso_lasso_coefs_result_fold_", fold_number, ".csv"))

  summary_test_ehr_model_coefficients <- read.csv(paste0("../models/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_ehr_model_result_fold_", fold_number, ".csv"))
  summary_test_null_model_coefficients <- read.csv(paste0("../models/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_null_model_result_fold_", fold_number, ".csv"))

  df <- readRDS("../data/t2d_cox_table.rds")
  covariates_df <- readRDS("../data/covariates.rds")
  df <- df %>% left_join(covariates_df)

  # Threshold events at 5 years
  df$Event <- sapply(1:nrow(df), function(i) {
    if (df$time_to_event[[i]] > 5) {
      0
    } else {
      df$Event[[i]]
    }
  })
  
  condition_df <- readRDS("../data/snomed_condition_frequency_table.rds")

  df$AGE <- df$age
  df$SEX <- df$sex

  model_means <- read.csv(paste0("../models/1k_snomed/1k_snomed_cox_lasso_model_means_result_fold_", fold_number, ".csv"),
                          row.names = "X")

  

  age_mean <- model_means["AGE", ]
  sex_mean <- model_means["SEX", ]
  prediction_mean <- model_means["prediction", ]

  basehazes <- read.csv(paste0("../models/1k_snomed/1k_snomed_cox_lasso_basehazes_result_fold_", fold_number, ".csv"))

  null_basehaz <- basehazes$null_basehaz
  ehr_basehaz <- basehazes$ehr_basehaz

  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(condition_df))

  ### SNOMED code imputation with nearest embedding neighbour
  all_embedding_terms <- row.names(embedding_table)

  lasso_coefs_absent <- lasso_coefs %>% filter(!(X %in% colnames(condition_df)))
  
  # absent_terms includes those in lasso_coefs not present in GS and present in embeddings_table
  absent_terms <- all_embedding_terms[all_embedding_terms %in% lasso_coefs_absent$X]
  absent_term_embeddings <- embedding_table[absent_terms, ]# lapply(lasso_coefs_absent$X, function(x) {get_embedding(x, embedding_table)})

  # Subset embedding table to only those present in target cohort
  target_cohort_embedding_terms <- all_embedding_terms[all_embedding_terms %in% colnames(condition_df)]
  target_cohort_embeddings <- embedding_table[target_cohort_embedding_terms, ] # lapply(gs_embedding_terms, function(x) {get_embedding(x, embedding_table)})

  # Find nearest neighbour in gs_embeddings for all absent lasso coef terms
  absent_term_nn_result <- nn2(data = target_cohort_embeddings, query = absent_term_embeddings, k = k)
  row.names(absent_term_nn_result$nn.idx) <- absent_terms

  # nn_terms is a matrix where rows correspond to lasso coef terms that are missing in GS. Number of columns corresponds to the number of nearest neighbours calculated.
  # Each entry is the jth nearest neighbour in GS for the ith missing lasso coef term, where i is the row number and j is the column number.
  nn_terms <- apply(absent_term_nn_result$nn.idx, c(1,2), function(x) {target_cohort_embedding_terms[x]})

  pre_imputation_df <- condition_df

  for (missing_term in row.names(nn_terms)) {
    missing_term_nns <- nn_terms[missing_term, ]
    nn_df <- condition_df[, missing_term_nns]
    nn_union <- as.numeric(apply(nn_df, 1, any))
    # Convert to sparse matrix so it can be cbinded with the condition_df
    nn_union <- as(as.matrix(nn_union), "sparseMatrix")
    colnames(nn_union) <- missing_term
    condition_df <- cbind(condition_df, nn_union)
  }

  # Recalculated lasso coefs present to check how many have been imputed
  pre_imputation_lasso_coefs_present <- lasso_coefs_present

  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(condition_df))
  ###

  age_centered <- df$age - age_mean
  sex_centered <- df$sex - sex_mean

  lasso_coef_names <- lasso_coefs_present$X
  ehr_coefs <- as.numeric(lasso_coefs_present$lasso_coefs)
  ehr_matrix <- as.matrix(condition_df[,lasso_coef_names])
  ehr_matrix <- +(ehr_matrix > 0)
  ehr_score <- ehr_matrix %*% ehr_coefs

  ehr_score_centered <- ehr_score - prediction_mean

  ehr_model_linear_predictor <- age_centered * (summary_test_ehr_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) + 
                                sex_centered * (summary_test_ehr_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                ehr_score_centered * (summary_test_ehr_model_coefficients %>% filter(X == "prediction") %>% pull(coef))
  ehr_survival_prediction <- exp(-ehr_basehaz) ^ exp(ehr_model_linear_predictor)
  ehr_onset_prediction <- 1 - ehr_survival_prediction

  ehr_roc <- pROC::roc(df$Event, ehr_onset_prediction)

  null_model_linear_predictor <- age_centered * (summary_test_null_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                              sex_centered * (summary_test_null_model_coefficients %>% filter(X == "SEX") %>% pull(coef))
  null_survival_prediction <- exp(-null_basehaz) ^ exp(null_model_linear_predictor)
  null_onset_prediction <- 1 - null_survival_prediction

  null_roc <- pROC::roc(df$Event, null_onset_prediction)

  missing_fg_snomed_codes <- setdiff(lasso_coefs$X, lasso_coefs_present$X)
  missing_lasso_coefs <- lasso_coefs[lasso_coefs$X %in% missing_fg_snomed_codes, ]

  rocs <- list(null = null_roc,
               ehr = ehr_roc)
  list(rocs = rocs,
       lasso_coefs = lasso_coefs, missing_lasso_coefs = missing_lasso_coefs)
}

results <- lapply(1:5, model_results_single_fold)

saveRDS(results, paste0("../results/1k_snomed/test_fg_model_fg_snomed_1k_", k, "nn_embedding_imputation_ukb_results.rds"))
