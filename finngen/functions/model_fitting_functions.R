library(glmnet)
library(dplyr)
library(survival)
library(Matrix)

train_and_test_survival_prs_interaction <- function(df, variables, training_index = NULL, alpha = 1, tte_threshold = NULL) {
  set.seed(42)
  
  if (is.null(training_index)) {
    
    training_set_size <- floor(4/5 * nrow(df))
    test_set_size <- nrow(df) - training_set_size
    training_index <- sample(c(rep(TRUE, training_set_size), rep(FALSE, test_set_size)))
    
  }
  training_df <- df[training_index, ]
  test_df <- df[!training_index, ]

  training_xs <- training_df[, variables]
  test_xs <- test_df[, variables]
  
  training_y <- pull(training_df, "T2D_INCLAVO")
  test_y <- pull(test_df, "T2D_INCLAVO")
  
  training_surv <- data.frame(time = training_df$time_to_event, status = training_y)
  test_surv <- data.frame(time = test_df$time_to_event, status = test_y)
  
  # Assume NAs come from joining condition and drug tables and can be treated as 0s
  training_xs <- training_xs %>% replace(is.na(.), 0)
  test_xs <- test_xs %>% replace(is.na(.), 0)
  
  gc()
  
  training_df <- NULL
  gc()
  

  glmnet_fold_ids <- get_stratified_k_fold_cv_indexes(training_surv,
                                                      k = 5,
                                                      id_colname = NULL,
                                                      event_colname = "status")
  
  lasso_model <- cv.glmnet(x = as.matrix(training_xs), y = as.matrix(training_surv), nfolds = 5, family = 'cox', trace.it = 1, foldid = glmnet_fold_ids, alpha = alpha) # , penalty.factor = penalty_factor)# , parallel = TRUE)
  
  
  lasso_coefs <- coef(lasso_model, s = 'lambda.min')
  lasso_coefs <- lasso_coefs[lasso_coefs[, 1] != 0, ]
  abs_lasso_coefs <- abs(lasso_coefs)
  lasso_coefs_sort_result <- sort(abs_lasso_coefs, decreasing = TRUE, index.return = TRUE)
  
  lasso_coefs <- lasso_coefs[lasso_coefs_sort_result$ix]
  
  lasso_variables <- names(lasso_coefs)
  
  training_xs_std_devs <- apply(training_xs, 2, sd)
  lasso_std_devs <- training_xs_std_devs[lasso_variables]
  scaled_lasso_coefs <- lasso_coefs * lasso_std_devs
  abs_scaled_lasso_coefs <- abs(scaled_lasso_coefs)
  scaled_lasso_coefs_sort_result <- sort(abs_scaled_lasso_coefs, decreasing = TRUE, index.return = TRUE)
  
  scaled_lasso_coefs <- scaled_lasso_coefs[scaled_lasso_coefs_sort_result$ix]
  
  n_cases_training <- sum(training_y == 1)
  
  model <- lasso_model
  
  training_prediction <- predict(model, as.matrix(training_xs), s = 'lambda.min', type = 'response')
  
  test_prediction <- predict(model, as.matrix(test_xs), s = 'lambda.min', type = 'response')
  
  test_df <- test_df %>% mutate(prediction = test_prediction)
  test_null_model <- coxph(Surv(time_to_event, T2D_INCLAVO) ~ AGE + SEX, data = test_df, method = "breslow")
  test_ehr_model <- coxph(Surv(time_to_event, T2D_INCLAVO) ~ AGE + SEX + prediction, data = test_df, method = "breslow")
  test_prs_model <- coxph(Surv(time_to_event, T2D_INCLAVO) ~ AGE + SEX + scale(SCORE1_AVG), data = test_df, method = "breslow")
  test_full_model <- coxph(Surv(time_to_event, T2D_INCLAVO) ~ AGE + SEX + prediction + scale(SCORE1_AVG), data = test_df, method = "breslow")
  test_interaction_model <- coxph(Surv(time_to_event, T2D_INCLAVO) ~ AGE + SEX + prediction*scale(SCORE1_AVG), data = test_df, method = "breslow")
  
  summary_test_null_model <- summary(test_null_model)
  summary_test_ehr_model <- summary(test_ehr_model)
  summary_test_prs_model <- summary(test_prs_model)
  summary_test_full_model <- summary(test_full_model)
  summary_test_interaction_model <- summary(test_interaction_model)
  
  null_survfit <- survfit(test_null_model)
  ehr_survfit <- survfit(test_ehr_model)
  prs_survfit <- survfit(test_prs_model)
  full_survfit <- survfit(test_full_model)
  interaction_survfit <- survfit(test_interaction_model)
  
  # each survfit$cumhaz has length equal to the number of distinct times. Index into the one corresponding to t == 5.
  null_basehaz <- null_survfit$cumhaz[[which(null_survfit$time >= 5)[[1]]]] 
  ehr_basehaz <- ehr_survfit$cumhaz[[which(ehr_survfit$time >= 5)[[1]]]] 
  prs_basehaz <- prs_survfit$cumhaz[[which(prs_survfit$time >= 5)[[1]]]] 
  full_basehaz <- full_survfit$cumhaz[[which(full_survfit$time >= 5)[[1]]]] 
  interaction_basehaz <- interaction_survfit$cumhaz[[which(interaction_survfit$time >= 5)[[1]]]] 
  
  null_prediction <- 1 - (exp(-null_basehaz) ^ exp(predict(test_null_model)))
  ehr_prediction <- 1 - (exp(-ehr_basehaz) ^ exp(predict(test_ehr_model)))
  prs_prediction <- 1 - (exp(-prs_basehaz) ^ exp(predict(test_prs_model)))
  full_prediction <- 1 - (exp(-full_basehaz) ^ exp(predict(test_full_model)))
  interaction_prediction <- 1 - (exp(-interaction_basehaz) ^ exp(predict(test_interaction_model)))
  
  # Threshold cases with TTE > tte_threshold
  if (!is.null(tte_threshold)) {
    test_df$T2D_INCLAVO <- ifelse(test_df$time_to_event > 5, 0, test_df$T2D_INCLAVO)
    test_y <- test_df$T2D_INCLAVO
  }

  null_roc <- pROC::roc(test_y, null_prediction)
  ehr_roc <- pROC::roc(test_y, ehr_prediction)
  prs_roc <- pROC::roc(test_y, prs_prediction)
  full_roc <- pROC::roc(test_y, full_prediction)
  interaction_roc <- pROC::roc(test_y, interaction_prediction)
  
  null_auc <- pROC::auc(test_y, null_prediction)
  ehr_auc <- pROC::auc(test_y, ehr_prediction)
  prs_auc <- pROC::auc(test_y, prs_prediction)
  full_auc <- pROC::auc(test_y, full_prediction)
  interaction_auc <- pROC::auc(test_y, interaction_prediction)
  
  summaries <- list(summary_test_null_model = summary_test_null_model,
                    summary_test_ehr_model = summary_test_ehr_model,
                    summary_test_prs_model = summary_test_prs_model,
                    summary_test_full_model = summary_test_full_model,
                    summary_test_interaction_model = summary_test_interaction_model)
  
  predictions <- list(null_prediction = null_prediction,
                      ehr_prediction = ehr_prediction,
                      prs_prediction = prs_prediction,
                      full_prediction = full_prediction,
                      interaction_prediction = interaction_prediction)
  
  rocs <- list(null_roc = null_roc,
               ehr_roc = ehr_roc,
               prs_roc = prs_roc,
               full_roc = full_roc,
               interaction_roc = interaction_roc)
  
  aucs <- list(null_auc = null_auc,
               ehr_auc = ehr_auc,
               prs_auc = prs_auc,
               full_auc = full_auc,
               interaction_auc = interaction_auc)
  
  basehazes <- list(null_basehaz = null_basehaz,
                    ehr_basehaz = ehr_basehaz,
                    prs_basehaz = ehr_basehaz,
                    full_basehaz = full_basehaz,
                    interaction_basehaz = interaction_basehaz)
  
  models <- list(null_model = test_null_model,
                 ehr_model = test_ehr_model,
                 prs_model = test_prs_model,
                 full_model = test_full_model,
                 interaction_model = test_interaction_model)
  
  list(model = model,
       lasso_coefs = lasso_coefs,
       scaled_lasso_coefs = scaled_lasso_coefs,
       training_surv = training_surv,
       test_surv = test_surv,
       training_index = training_index,
       training_prediction = training_prediction,
       test_prediction = test_prediction,
       summaries = summaries,
       predictions = predictions,
       rocs = rocs,
       aucs = aucs,
       basehazes = basehazes,
       cox_models = models
  )
}


get_k_fold_cv_indexes <- function(df, k, sample_seed = 42) {
  number_of_rows <- nrow(df)
  cv_indexes <- rep(1:k, length.out = number_of_rows)
  set.seed(sample_seed)
  sample(cv_indexes)
}

# Need to remove id_colname as it is not used
get_stratified_k_fold_cv_indexes <- function(df, k, id_colname, event_colname, sample_seed = 42) {
  # Create temporary index column to reorder the final cv indexes
  df$tmp_original_index <- 1:nrow(df)
  cases <- df[df[, event_colname] == 1, ]
  controls <- df[df[, event_colname] == 0, ]
  df <- rbind(cases, controls)
  case_indexes <- get_k_fold_cv_indexes(cases, k = k, sample_seed = sample_seed)
  control_indexes <- get_k_fold_cv_indexes(controls, k = k, sample_seed = sample_seed + 1)
  cv_indexes <- c(case_indexes, control_indexes)
  
  # Calculate the index required to reorder cv_indexes according to the original df order
  cv_indexes_reorder_index <- order(df$tmp_original_index, decreasing = FALSE)
  cv_indexes <- cv_indexes[cv_indexes_reorder_index]
}


get_df <- function() {
  condition_preprocessing_result <- readRDS("data/r11_condition/20231112/condition_preprocessing_result.rds")
  
  t2d_endpoint <- condition_preprocessing_result$t2d_endpoint
  
  condition_id_name_table <- condition_preprocessing_result$condition_id_name_table
  
  condition_presence_wide <- condition_preprocessing_result$condition_presence_wide
  
  condition_preprocessing_result <- NULL
  gc()
  
  variables <- condition_presence_wide %>% select(-person_id) %>% colnames
  
  df <- t2d_endpoint %>% inner_join(condition_presence_wide)
  
  condition_presence_wide <- NULL
  t2d_endpoint <- NULL
  gc()
  list(df = df, penalised_variables = variables)
}


get_abs_sorted_coefs <- function(model) {
  coefs <- coef(model, s = 'lambda.1se')
  coefs <- coefs[coefs[, 1] != 0, ]
  
  abs_coefs <- abs(coefs)
  coefs_sort_result <- sort(abs_coefs, decreasing = TRUE, index.return = TRUE)
  
  coefs[coefs_sort_result$ix]
}