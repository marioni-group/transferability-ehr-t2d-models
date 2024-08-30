library(dplyr)
library(glmnet)

export_lasso_and_incremental_coefs <- function(result_file, result_path, output_path, model_description) {

  result <- readRDS(paste0(result_path, result_file, ".rds"))
  
  write_csv_with_result_suffix <- function(x, prefix) {
    write.csv(x, paste0(output_path, model_description, "_", prefix, "_", result_file, ".csv"))
  }
  
  incremental_coefficients <- lapply(result$summaries, function(x){x$coefficients})
  
  lasso_coefs <- result$lasso_coefs
  
  coef_list <- list(lasso = lasso_coefs, incremental = incremental_coefficients)
  
  sapply(names(incremental_coefficients), function(x) {
    write_csv_with_result_suffix(incremental_coefficients[[x]], paste0("incremental_coefs_", x))
  })
  
  write_csv_with_result_suffix(as.data.frame(lasso_coefs), "lasso_coefs")
  
  
  interaction_model_means <- result$cox_models$interaction_model$means
  write_csv_with_result_suffix(as.data.frame(interaction_model_means), "model_means")
  
  basehazes <- result$basehazes
  write_csv_with_result_suffix(as.data.frame(basehazes), "basehazes")
}

for (cl_i in 1:5) {
    export_lasso_and_incremental_coefs(result_file = paste0("result_cl_fold_", cl_i, "_tm_fold_1"),
                                       result_path = "models/r11_condition/20240814_topic_model_all_folds/cox_lasso_models/",
                                       output_path = "models/r11_condition/20240814_topic_model_all_folds/cox_lasso_models/extracted_model_tables/",
                                       model_description = "snomed_grouped_by_topic")
}

  post_processing_result <- readRDS(paste0("models/r11_condition/20240814_topic_model_all_folds/process_topic_model_result_200_topics_all_folds.rds"))
  condition_topics_unique <- post_processing_result$condition_topics_unique
  write.csv(condition_topics_unique, paste0("models/r11_condition/20240814_topic_model_all_folds/cox_lasso_models/extracted_model_tables/condition_topics_unique.csv"))
