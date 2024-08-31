library(glmnet)
library(dplyr)
# library(SurvMetrics)
library(survival)

model_files_folder <- "/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/topic_model_grouping_tm_all_folds/"

model_results_single_fold <- function(tm_fold_number, cl_fold_number) {

  lasso_coefs <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_lasso_coefs_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_ehr_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_ehr_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_full_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_full_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_interaction_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_interaction_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_null_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_null_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  summary_test_prs_model_coefficients <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_incremental_coefs_summary_test_prs_model_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  
  df <- readRDS(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table_snomed_grouped_by_topic_model_tm_all_folds_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".rds"))
  
  # Remove participants with TTE == 0
  df <- df %>% filter(time_to_event > 0)


  t2d_prs <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/INTERVENE/flagship/prs_output/T2D_PRS.sscore", sep = "\t")
  t2d_prs$IID <- as.character(t2d_prs$IID)
  # Add PRS
  df <- df %>% left_join(t2d_prs, by = c("id" = "IID"))
  df <- df %>% filter(!is.na(SCORE1_SUM))
  # Scale PRS
  df$SCORE1_SUM <- scale(df$SCORE1_SUM)
  
  means <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_model_means_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  
  age_mean <- means %>% filter(X == "AGE") %>% pull(interaction_model_means)
  sex_mean <- means %>% filter(X == "SEX") %>% pull(interaction_model_means)
  prediction_mean <- means %>% filter(X == "prediction") %>% pull(interaction_model_means)
  scaled_prs_mean <- means %>% filter(X == "scale(SCORE1_AVG)") %>% pull(interaction_model_means)
  prediction_scaled_prs_interaction_mean <- means %>% filter(X == "prediction:scale(SCORE1_AVG)") %>% pull(interaction_model_means)
  
  basehazes <- read.csv(paste0(model_files_folder, "snomed_grouped_by_topic_basehazes_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))
  
  null_basehaz <- basehazes$null_basehaz
  ehr_basehaz <- basehazes$ehr_basehaz
  prs_basehaz <- basehazes$prs_basehaz
  full_basehaz <- basehazes$full_basehaz
  interaction_basehaz <- basehazes$interaction_basehaz
  
  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(df))
  
  df$age_centered <- df$Age - age_mean
  df$sex_centered <- df$sex - sex_mean
  
  lasso_coef_names <- lasso_coefs_present$X
  ehr_coefs <- as.numeric(lasso_coefs_present$lasso_coefs)
  ehr_matrix <- as.matrix(df[,lasso_coef_names])
  ehr_score <- ehr_matrix %*% ehr_coefs
  
  df$ehr_score_centered <- ehr_score - prediction_mean
  
  df$scaled_prs_centered <- df$SCORE1_SUM - scaled_prs_mean
  
  df$interaction_centered <- ehr_score * df$SCORE1_SUM - prediction_scaled_prs_interaction_mean
  
  if (filter_by_quantile) {
    # Filter by variable quantile
    df <- df %>% mutate(var_quantile = ntile(!!ensym(quantile_variable), n_quantiles))
    df <- df %>% filter(var_quantile == filter_quantile)
  }
  
  if (filter_cases_gt_2_years) {
    df <- df %>% filter((Event == 0) | (time_to_event > 2))
  }
  
  ehr_model_linear_predictor <- df$age_centered * (summary_test_ehr_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) + 
                                df$sex_centered * (summary_test_ehr_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                df$ehr_score_centered * (summary_test_ehr_model_coefficients %>% filter(X == "prediction") %>% pull(coef))
  ehr_survival_prediction <- exp(-ehr_basehaz) ^ exp(ehr_model_linear_predictor)
  ehr_onset_prediction <- 1 - ehr_survival_prediction
  
  ehr_roc <- pROC::roc(df$Event, ehr_onset_prediction)
  ehr_prauc <- MLmetrics::PRAUC(ehr_onset_prediction, df$Event)
  
  null_model_linear_predictor <- df$age_centered * (summary_test_null_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                 df$sex_centered * (summary_test_null_model_coefficients %>% filter(X == "SEX") %>% pull(coef))
  null_survival_prediction <- exp(-null_basehaz) ^ exp(null_model_linear_predictor)
  null_onset_prediction <- 1 - null_survival_prediction
  
  null_roc <- pROC::roc(df$Event, null_onset_prediction)
  null_prauc <- MLmetrics::PRAUC(null_onset_prediction, df$Event)
  
  prs_model_linear_predictor <- df$age_centered * (summary_test_prs_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                df$sex_centered * (summary_test_prs_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                df$scaled_prs_centered * (summary_test_prs_model_coefficients %>% filter(X == "scale(SCORE1_AVG)") %>% pull(coef))
  prs_survival_prediction <- exp(-prs_basehaz) ^ exp(prs_model_linear_predictor)
  prs_onset_prediction <- 1 - prs_survival_prediction
  
  prs_roc <- pROC::roc(df$Event, prs_onset_prediction)
  prs_prauc <- MLmetrics::PRAUC(prs_onset_prediction, df$Event)
  
  full_model_linear_predictor <- df$age_centered * (summary_test_full_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                 df$sex_centered * (summary_test_full_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                 df$scaled_prs_centered * (summary_test_full_model_coefficients %>% filter(X == "scale(SCORE1_AVG)") %>% pull(coef)) +
                                 df$ehr_score_centered * (summary_test_full_model_coefficients %>% filter(X == "prediction") %>% pull(coef))
  full_survival_prediction <- exp(-full_basehaz) ^ exp(full_model_linear_predictor)
  full_onset_prediction <- 1 - full_survival_prediction
  
  full_roc <- pROC::roc(df$Event, full_onset_prediction)
  full_prauc <- MLmetrics::PRAUC(full_onset_prediction, df$Event)
  
  interaction_model_linear_predictor <- df$age_centered * (summary_test_interaction_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                        df$sex_centered * (summary_test_interaction_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                        df$scaled_prs_centered * (summary_test_interaction_model_coefficients %>% filter(X == "scale(SCORE1_AVG)") %>% pull(coef)) +
                                        df$ehr_score_centered * (summary_test_interaction_model_coefficients %>% filter(X == "prediction") %>% pull(coef)) +
                                        df$interaction_centered * (summary_test_interaction_model_coefficients %>% filter(X == "prediction:scale(SCORE1_AVG)") %>% pull(coef))
  interaction_survival_prediction <- exp(-interaction_basehaz) ^ exp(interaction_model_linear_predictor)
  interaction_onset_prediction <- 1 - interaction_survival_prediction
  
  interaction_roc <- pROC::roc(df$Event, interaction_onset_prediction)
  interaction_prauc <- MLmetrics::PRAUC(interaction_onset_prediction, df$Event)
 
  missing_fg_snomed_codes <- setdiff(lasso_coefs$X, lasso_coefs_present$X)
  missing_lasso_coefs <- lasso_coefs[lasso_coefs$X %in% missing_fg_snomed_codes, ]

  rocs <- list(null = null_roc,
               ehr = ehr_roc,
               prs = prs_roc,
               full = full_roc,
               interaction = interaction_roc)

  aucs <- sapply(rocs, function(x) {x$auc})


  brier_times <- 1:5

  list(aucs = aucs, rocs = rocs, 
       lasso_coefs = lasso_coefs, missing_lasso_ceofs = missing_lasso_coefs)
}

results <- lapply(1, function(tm_fold_i) {
  lapply(1:5, function(cl_fold_i) {
    model_results_single_fold(tm_fold_i, cl_fold_i)
  })
})

saveRDS(results, "/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/results/r11/test_fg_model_fg_snomed_grouped_by_topic_tm_all_folds_results_with_auc.rds")
