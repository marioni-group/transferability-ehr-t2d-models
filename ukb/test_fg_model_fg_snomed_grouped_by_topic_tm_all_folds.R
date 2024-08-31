library(glmnet)
library(dplyr)
library(survival)


model_files_folder <- "../models/topic_model_grouping_tm_all_folds/"

model_results_single_fold <- function(tm_fold_number, cl_fold_number) {

  lasso_coefs <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_lasso_coefs_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_ehr_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_ehr_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_null_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_null_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))

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
  
  condition_df <- readRDS(paste0("../data/presence_snomed_grouped_by_topic_model_tm_all_folds_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".rds"))
  
  df <- cbind(df, condition_df)

 
  means <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_model_means_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  
  age_mean <- means %>% filter(X == "AGE") %>% pull(interaction_model_means)
  sex_mean <- means %>% filter(X == "SEX") %>% pull(interaction_model_means)
  prediction_mean <- means %>% filter(X == "prediction") %>% pull(interaction_model_means)

  basehazes <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_basehazes_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  
  null_basehaz <- basehazes$null_basehaz
  ehr_basehaz <- basehazes$ehr_basehaz

  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(df))
  
  df$age_centered <- df$age - age_mean
  df$sex_centered <- df$sex - sex_mean
  
  lasso_coef_names <- lasso_coefs_present$X
  ehr_coefs <- as.numeric(lasso_coefs_present$lasso_coefs)
  ehr_matrix <- as.matrix(df[,lasso_coef_names])
  ehr_score <- ehr_matrix %*% ehr_coefs
  
  df$ehr_score_centered <- ehr_score - prediction_mean
  
  ehr_model_linear_predictor <- df$age_centered * (summary_test_ehr_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) + 
                                df$sex_centered * (summary_test_ehr_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                df$ehr_score_centered * (summary_test_ehr_model_coefficients %>% filter(X == "prediction") %>% pull(coef))
  ehr_survival_prediction <- exp(-ehr_basehaz) ^ exp(ehr_model_linear_predictor)
  ehr_onset_prediction <- 1 - ehr_survival_prediction
  
  ehr_roc <- pROC::roc(df$Event, ehr_onset_prediction)

  null_model_linear_predictor <- df$age_centered * (summary_test_null_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                 df$sex_centered * (summary_test_null_model_coefficients %>% filter(X == "SEX") %>% pull(coef))
  null_survival_prediction <- exp(-null_basehaz) ^ exp(null_model_linear_predictor)
  null_onset_prediction <- 1 - null_survival_prediction
  
  null_roc <- pROC::roc(df$Event, null_onset_prediction)
  
  missing_fg_snomed_codes <- setdiff(lasso_coefs$X, lasso_coefs_present$X)
  missing_lasso_coefs <- lasso_coefs[lasso_coefs$X %in% missing_fg_snomed_codes, ]
 
  
  rocs <- list(null = null_roc,
               ehr = ehr_roc) 
  
  list(rocs = rocs, lasso_coefs = lasso_coefs, missing_lasso_coefs = missing_lasso_coefs)
}

results <- lapply(1, function(tm_fold_i) {
  lapply(1:5, function(cl_fold_i) {
    model_results_single_fold(tm_fold_i, cl_fold_i)
  })
})

saveRDS(results, "../results/test_fg_model_fg_snomed_grouped_by_topic_tm_all_folds_results.rds")
