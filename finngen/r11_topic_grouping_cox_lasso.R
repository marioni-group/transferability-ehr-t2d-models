library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(glmnet)
library(randomForest)
library(ggpubr)
source('src/functions/prs_functions.R')
source('src/functions/model_fitting_functions.R')
source('src/functions/util_functions.R')

topic_grouping_cox_lasso_single_fold <- function(topic_model_fold_number, cox_lasso_fold_number) {

  data_folder_path <- "data/r11_condition/20231113/"
  topic_model_path <- "models/r11_condition/20231124/topic_models/"
  save_path <- "models/r11_condition/20231124/cox_lasso_models/"
  
  condition_preprocessing_result <- read_rds(paste0(data_folder_path, "condition_preprocessing_result.rds"))
  
  t2d_endpoint <- condition_preprocessing_result$t2d_endpoint
  
  condition_presence_wide <- condition_preprocessing_result$condition_presence_wide
  
  condition_preprocessing_result <- NULL
  gc()
  
  t2d_endpoint <- t2d_endpoint %>% filter(prediction_model_set_index == 1)
  condition_presence_wide <- condition_presence_wide %>% filter(person_id %in% t2d_endpoint$person_id)
  (t2d_endpoint$person_id == condition_presence_wide$person_id) %>% all

  gc()
  
  variables <- colnames(condition_presence_wide %>% select(-person_id))
  
  # Import term-topic table and use to group terms
  condition_topics_unique <- readRDS(paste0(topic_model_path, "min_ppxt_post_processing_result_fold_", topic_model_fold_number, ".rds"))$condition_topics_unique
  topic_presence_table <- matrix(data = 0, nrow = nrow(condition_presence_wide), ncol = max(condition_topics_unique$topic))

  # Loop through topic numbers. 
  # For each topic, for each term in the topic:
  # Take the OR of the corresponding column in condition_presence_wide and the topic column in topic_presence_table
  for (topic_number in 1:ncol(topic_presence_table)) {
    terms <- condition_topics_unique %>% filter(topic == topic_number) %>% pull(term)
    for (term in terms) {
      if (term %in% colnames(condition_presence_wide)) {
        topic_presence_table[, topic_number ] <- topic_presence_table[, topic_number] | condition_presence_wide[, term]
      }
    }
  }
  
  colnames(topic_presence_table) <- paste0("t", 1:ncol(topic_presence_table))
  
  stop()
  
  variables <- colnames(topic_presence_table)
  
  condition_presence_table <- bind_cols(list(t2d_endpoint, topic_presence_table))
  
  # Add PRS
  condition_presence_table <- add_t2d_scores_to_df(condition_presence_table, 'FINNGENID')$df
  
  
  k_fold_cv_indexes <- t2d_endpoint$prediction_model_set_k_fold_cv_indexes
  training_index <- !(k_fold_cv_indexes == cox_lasso_fold_number)
  result <- train_and_test_survival_prs_interaction(condition_presence_table, variables, training_index = training_index)
  
  lasso_coef_names <- names(result$lasso_coefs)
  
  topic_model_post_processing_result <- readRDS(paste0(topic_model_path, "min_ppxt_post_processing_result_fold_", topic_model_fold_number, ".rds"))
  
  plot_selected_topics_top_terms <- function(post_processing_result, ts, save_folder_path = results_folder_path, lasso_betas = result$lasso_coefs) {
    lapply(ts, function(t) {
      post_processing_result$ap_top_terms %>%
        filter(topic == t) %>%
        mutate(term = reorder_within(term, beta, topic)) %>%
        ggplot(aes(beta, term)) +
        geom_col(show.legend = FALSE, width = 0.2) +
        ggtitle(paste0("Topic ", t, ", lasso coef = ", round(lasso_betas[[paste0("t", t)]], digits = 2))) +
        scale_y_reordered() +
        theme_minimal(base_size = 16) + theme(aspect.ratio = 1/2)
    })
  }
  
  top_topics_top_term_plots <- plot_selected_topics_top_terms(topic_model_post_processing_result,
                                                              ts = sapply(lasso_coef_names[1:9], function(x) {substr(x, 2, nchar(x))}))
  
  ggarrange(plotlist = top_topics_top_term_plots, ncol = 3, nrow = 3)
  
  saveRDS(result, paste0(save_path, "result_cl_fold_", cox_lasso_fold_number, "_tm_fold_", topic_model_fold_number, ".rds"))
}

for (cl_i in 1:5) {
  for (tm_i in 1:3) {
    print(paste0("Current cox lasso fold: ", cl_i))
    print(paste0("Current topic model fold: ", tm_i))
    system.time({topic_grouping_cox_lasso_single_fold(topic_model_fold_number = tm_i,
                                         cox_lasso_fold_number = cl_i)})
    gc()
  }
}

cl_i <- NULL
tm_i <- NULL

aucs <- lapply(1:5, function(cl_i) {
  sapply(1:3, function(tm_i) {
    result <- readRDS(paste0("models/r11_condition/20231124/cox_lasso_models/result_cl_fold_", cl_i, "_tm_fold_", tm_i, ".rds"))
    result$aucs
  })
})