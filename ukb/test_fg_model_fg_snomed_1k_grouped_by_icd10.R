library(glmnet)
library(dplyr)
library(survival)

model_results_single_fold <- function(fold_number) {
  model_files_path <- "../models/icd_category_grouping/"

  lasso_coefs <- read.csv(paste0(model_files_path, "snomed_grouped_by_icd_1k_freq_cox_lasso_lasso_coefs_result_fold_", fold_number, ".csv"))

  summary_test_ehr_model_coefficients <- read.csv(paste0("../models/icd_category_grouping/snomed_grouped_by_icd_1k_freq_cox_lasso_incremental_coefs_summary_test_ehr_model_result_fold_", fold_number, ".csv"))
  summary_test_null_model_coefficients <- read.csv(paste0("../models/icd_category_grouping/snomed_grouped_by_icd_1k_freq_cox_lasso_incremental_coefs_summary_test_null_model_result_fold_", fold_number, ".csv"))

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

  condition_df <- readRDS(paste0("../data/presence_snomed_grouped_by_icd10_fold_", fold_number, ".rds"))

  df <- cbind(df, condition_df)

  df$AGE <- df$age
  df$SEX <- df$sex

  model_means <- read.csv(paste0("../models/icd_category_grouping/snomed_grouped_by_icd_1k_freq_cox_lasso_model_means_result_fold_", fold_number, ".csv"),
                          row.names = "X")

  

  age_mean <- model_means["AGE", ]
  sex_mean <- model_means["SEX", ]
  prediction_mean <- model_means["prediction", ]

  basehazes <- read.csv(paste0("../models/icd_category_grouping/snomed_grouped_by_icd_1k_freq_cox_lasso_basehazes_result_fold_", fold_number, ".csv"))

  null_basehaz <- basehazes$null_basehaz
  ehr_basehaz <- basehazes$ehr_basehaz

  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(df))

  age_centered <- df$age - age_mean
  sex_centered <- df$sex - sex_mean

  lasso_coef_names <- lasso_coefs_present$X
  ehr_coefs <- as.numeric(lasso_coefs_present$lasso_coefs)
  ehr_matrix <- as.matrix(df[,lasso_coef_names])
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

  null_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = null_survival_prediction) # , t_star = 5)
  ehr_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = ehr_survival_prediction) # , t_star = 5)

  null_brier_score <- Brier(Surv(df$time_to_event), predicted = null_survival_prediction, t_star = 5)
  ehr_brier_score <- Brier(Surv(df$time_to_event), predicted = ehr_survival_prediction, t_star = 5)

  rocs <- list(null = null_roc,
               ehr = ehr_roc)

  aucs <- sapply(rocs, function(x) {
    x$auc
  })

  c_indexes <- list(null = null_c_index,
                    ehr = ehr_c_index)

  brier_scores <- list(null = null_brier_score,
                       ehr = ehr_brier_score)


  list(aucs = aucs,
       c_indexes = c_indexes,
       brier_scores = brier_scores,
       lasso_coefs = lasso_coefs, missing_lasso_coefs = missing_lasso_coefs)
}

results <- lapply(1:5, model_results_single_fold)

saveRDS(results, "../results/test_fg_model_fg_snomed_1k_grouped_by_icd10_auc_cindex_brier.rds")
