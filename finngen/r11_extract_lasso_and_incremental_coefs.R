library(dplyr)
library(glmnet)

result_file <- "result_fold_5"
result <- readRDS(paste0("models/r11_condition/20240223/1k_snomed_cox_lasso/", result_file, ".rds"))

output_path <- "models/r11_condition/20240223/1k_snomed_cox_lasso/extracted_model_tables/"
model_description <- "1k_snomed_cox_lasso"

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
