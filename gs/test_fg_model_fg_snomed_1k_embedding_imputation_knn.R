library(glmnet)
library(dplyr)
library(SurvMetrics)
library(survival)
library(RANN)
source("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/src/r11/snomed_embedding_functions.R")
embedding_table <- load_embedding_table()

model_results_single_fold <- function(fold_number, k) {
  model_files_path <- "/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/"

  lasso_coefs <- read.csv(paste0(model_files_path, "1k_snomed_cox_lasso_lasso_coefs_result_fold_", fold_number, ".csv"))

  summary_test_ehr_model_coefficients <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_ehr_model_result_fold_", fold_number, ".csv"))
  summary_test_full_model_coefficients <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_full_model_result_fold_", fold_number, ".csv"))
  summary_test_interaction_model_coefficients <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_interaction_model_result_fold_", fold_number, ".csv"))
  summary_test_null_model_coefficients <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_null_model_result_fold_", fold_number, ".csv"))
  summary_test_prs_model_coefficients <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_incremental_coefs_summary_test_prs_model_result_fold_", fold_number, ".csv"))
  
  df <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table.rds")
  
    # Remove participants with TTE == 0
  df <- df %>% filter(time_to_event > 0)
  

  df$AGE <- df$Age
  df$SEX <- df$sex
  t2d_prs <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/INTERVENE/flagship/prs_output/T2D_PRS.sscore", sep = "\t")
  t2d_prs$IID <- as.character(t2d_prs$IID)
  # Add PRS
  df <- df %>% left_join(t2d_prs, by = c("id" = "IID"))
  df <- df %>% filter(!is.na(SCORE1_SUM))
  # Scale PRS
  df$SCORE1_SUM <- scale(df$SCORE1_SUM)

  model_means <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_model_means_result_fold_", fold_number, ".csv"),
                          row.names = "X")

  

  age_mean <- model_means["AGE", ]
  sex_mean <- model_means["SEX", ]
  prediction_mean <- model_means["prediction", ]
  scaled_prs_mean <- model_means["scale(SCORE1_AVG)", ]
  prediction_scaled_prs_interaction_mean <- model_means["prediction:scale(SCORE1_AVG)", ]

  basehazes <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_basehazes_result_fold_", fold_number, ".csv"))

  null_basehaz <- basehazes$null_basehaz
  ehr_basehaz <- basehazes$ehr_basehaz
  prs_basehaz <- basehazes$prs_basehaz
  full_basehaz <- basehazes$full_basehaz
  interaction_basehaz <- basehazes$interaction_basehaz

  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(df))

  # SNOMED code imputation with nearest embedding neighbour
  all_embedding_terms <- row.names(embedding_table)

  lasso_coefs_absent <- lasso_coefs %>% filter(!(X %in% colnames(df)))
  
  # absent_terms includes those in lasso_coefs not present in GS and present in embeddings_table
  absent_terms <- all_embedding_terms[all_embedding_terms %in% lasso_coefs_absent$X]
  absent_term_embeddings <- embedding_table[absent_terms, ]# lapply(lasso_coefs_absent$X, function(x) {get_embedding(x, embedding_table)})
  # absent_term_embedding_df <- do.call(rbind, absent_term_embeddings)

  # Subset embedding table to only those present in GS
  gs_embedding_terms <- all_embedding_terms[all_embedding_terms %in% colnames(df)]
  gs_embeddings <- embedding_table[gs_embedding_terms, ] # lapply(gs_embedding_terms, function(x) {get_embedding(x, embedding_table)}) 

  # Find nearest neighbour in gs_embeddings for all absent lasso coef terms
  absent_term_nn_result <- nn2(data = gs_embeddings, query = absent_term_embeddings, k = k)
  row.names(absent_term_nn_result$nn.idx) <- absent_terms


  # nn_terms is a matrix where rows correspond to lasso coef terms that are missing in GS. Number of columns corresponds to the number of nearest neighbours calculated.
  # Each entry is the jth nearest neighbour in GS for the ith missing lasso coef term, where i is the row number and j is the column number.
  nn_terms <- apply(absent_term_nn_result$nn.idx, c(1,2), function(x) {gs_embedding_terms[x]})


  # Impute missing cols in df with the nearest neigbour

  pre_imputation_df <- df

  for (missing_term in row.names(nn_terms)) {
    missing_term_nns <- nn_terms[missing_term, ]
    nn_df <- df[, missing_term_nns]
    nn_union <- as.numeric(apply(nn_df, 1, any))
    df[[missing_term]] <- nn_union
  }
  

  # Recalculated lasso coefs present to check how many have been imputed
  pre_imputation_lasso_coefs_present <- lasso_coefs_present

  lasso_coefs_present <- lasso_coefs %>% filter(X %in% colnames(df))

  age_centered <- df$Age - age_mean
  sex_centered <- df$sex - sex_mean

  lasso_coef_names <- lasso_coefs_present$X
  ehr_coefs <- as.numeric(lasso_coefs_present$lasso_coefs)
  ehr_matrix <- as.matrix(df[,lasso_coef_names])
  ehr_score <- ehr_matrix %*% ehr_coefs

  ehr_score_centered <- ehr_score - prediction_mean
  scaled_prs_centered <- df$SCORE1_SUM - scaled_prs_mean

  interaction_centered <- ehr_score * df$SCORE1_SUM - prediction_scaled_prs_interaction_mean

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

  prs_model_linear_predictor <- age_centered * (summary_test_prs_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                sex_centered * (summary_test_prs_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                scaled_prs_centered * (summary_test_prs_model_coefficients %>% filter(X == "scale(SCORE1_AVG)") %>% pull(coef))
  prs_survival_prediction <- exp(-prs_basehaz) ^ exp(prs_model_linear_predictor)
  prs_onset_prediction <- 1 - prs_survival_prediction

  prs_roc <- pROC::roc(df$Event, prs_onset_prediction)

  full_model_linear_predictor <- age_centered * (summary_test_full_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                 sex_centered * (summary_test_full_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                 scaled_prs_centered * (summary_test_full_model_coefficients %>% filter(X == "scale(SCORE1_AVG)") %>% pull(coef)) +
                                 ehr_score_centered * (summary_test_full_model_coefficients %>% filter(X == "prediction") %>% pull(coef))
  full_survival_prediction <- exp(-full_basehaz) ^ exp(full_model_linear_predictor)
  full_onset_prediction <- 1 - full_survival_prediction

  full_roc <- pROC::roc(df$Event, full_onset_prediction)

  interaction_model_linear_predictor <- age_centered * (summary_test_interaction_model_coefficients %>% filter(X == "AGE") %>% pull(coef)) +
                                        sex_centered * (summary_test_interaction_model_coefficients %>% filter(X == "SEX") %>% pull(coef)) +
                                        scaled_prs_centered * (summary_test_interaction_model_coefficients %>% filter(X == "scale(SCORE1_AVG)") %>% pull(coef)) +
                                        ehr_score_centered * (summary_test_interaction_model_coefficients %>% filter(X == "prediction") %>% pull(coef)) +
                                        interaction_centered * (summary_test_interaction_model_coefficients %>% filter(X == "prediction:scale(SCORE1_AVG)") %>% pull(coef))
  interaction_survival_prediction <- exp(-interaction_basehaz) ^ exp(interaction_model_linear_predictor)
  interaction_onset_prediction <- 1 - interaction_survival_prediction

  interaction_roc <- pROC::roc(df$Event, interaction_onset_prediction)

  missing_fg_snomed_codes <- setdiff(lasso_coefs$X, lasso_coefs_present$X)
  missing_lasso_coefs <- lasso_coefs[lasso_coefs$X %in% missing_fg_snomed_codes, ]

  null_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = null_survival_prediction) # , t_star = 5)
  ehr_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = ehr_survival_prediction) # , t_star = 5)
  prs_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = prs_survival_prediction) # , t_star = 5)
  full_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = full_survival_prediction) # , t_star = 5)
  interaction_c_index <- Cindex(Surv(df$time_to_event, df$Event), predicted = interaction_survival_prediction) # , t_star = 5)

  rocs <- list(null = null_roc,
               ehr = ehr_roc,
               prs = prs_roc,
               full = full_roc,
               interaction = interaction_roc)

  aucs <- sapply(rocs, function(x) {
    x$auc
  })

  c_indexes <- list(null = null_c_index,
                    ehr = ehr_c_index,
                    prs = prs_c_index,
                    full = full_c_index,
                    interaction = interaction_c_index)

  brier_times <- 1:5

  null_briers <- sapply(brier_times, function(b_t) {
    Brier(Surv(df$time_to_event, df$Event), pre_sp = null_survival_prediction, t_star = b_t)
  })
  ehr_briers <- sapply(brier_times, function(b_t) {
    Brier(Surv(df$time_to_event, df$Event), pre_sp = ehr_survival_prediction, t_star = b_t)
  })
  prs_briers <- sapply(brier_times, function(b_t) {
    Brier(Surv(df$time_to_event, df$Event), pre_sp = prs_survival_prediction, t_star = b_t)
  })
  full_briers <- sapply(brier_times, function(b_t) {
    Brier(Surv(df$time_to_event, df$Event), pre_sp = full_survival_prediction, t_star = b_t)
  })
  interaction_briers <- sapply(brier_times, function(b_t) {
    Brier(Surv(df$time_to_event, df$Event), pre_sp = interaction_survival_prediction, t_star = b_t)
  })

  briers_list <- list(
    null = null_briers,
    ehr = ehr_briers,
    prs = prs_briers,
    full = full_briers,
    interaction = interaction_briers
  )

  briers <- as.data.frame(briers_list)

  list(aucs = aucs, c_indexes = c_indexes, briers = briers,
       lasso_coefs = lasso_coefs, missing_lasso_ceofs = missing_lasso_coefs)
}

for (k_value in c(10, 20, 50, 100)) {
  print(system.time({
    results <- lapply(1:5, model_results_single_fold, k = k_value)
  }))
  saveRDS(results, paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/results/r11/test_fg_model_fg_snomed_1k_embedding_imputation_", k_value, "nn_results_with_auc_cindex_brier.rds"))
}