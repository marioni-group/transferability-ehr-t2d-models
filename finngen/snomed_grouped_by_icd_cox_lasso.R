library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(glmnet)
library(randomForest)
source('src/functions/prs_functions.R')
source('src/functions/model_fitting_functions.R')
source('src/functions/util_functions.R')
source('src/functions/r11_group_snomed_by_icd.R')


data_path <- "data/r11_condition/20231113/"
save_path <- "models/r11_condition/20231115/snomed_grouped_by_icd_1k_freq_cox_lasso/"
fold_number <- 1

load_file <- function(filename) {
  if (endsWith(tolower(filename), ".rds")) {
    readRDS(paste0(data_path, filename))
  } else {
    stop(paste0("Error: unsupported file extension in load_file: ", filename))
  }
}

condition_preprocessing_result <- load_file("condition_preprocessing_result.rds")

t2d_endpoint <- condition_preprocessing_result$t2d_endpoint
gc()

condition_id_name_table <- condition_preprocessing_result$condition_id_name_table

condition_presence_wide <- condition_preprocessing_result$condition_presence_wide

condition_preprocessing_result <- NULL
gc()

# Filter files to prediction model set
t2d_endpoint <- t2d_endpoint %>% filter(prediction_model_set_index == 1)
condition_presence_wide <- condition_presence_wide %>% filter(person_id %in% t2d_endpoint$person_id)
(t2d_endpoint$person_id == condition_presence_wide$person_id) %>% all
gc()

k_fold_cv_indexes <- t2d_endpoint$prediction_model_set_k_fold_cv_indexes

training_condition_presence_wide <- condition_presence_wide[!(k_fold_cv_indexes == fold_number), ]

condition_counts <- training_condition_presence_wide %>% select(-person_id) %>% apply(2, sum) %>% sort(decreasing = TRUE)

training_condition_presence_wide <- NULL

gc()

conditions <- condition_counts %>% names

condition_presence_wide <- condition_presence_wide %>% select(all_of(conditions))
gc()

variables <- conditions %>% sample


################################################################################
# Grouping by ICD10 code

icd_presence_wide <- setNames(data.frame(matrix(0, ncol = length(unique_mapped_icd10), nrow = nrow(condition_presence_wide))), unique_mapped_icd10)

for (variable in variables) {
  icd_codes <- ifelse(variable %in% names(snomed_to_icd10_lookup), snomed_to_icd10_lookup[[variable]], "")
  icd_codes <- strsplit(icd_codes, split = ",")[[1]]
  if (length(icd_codes) > 0) {
    for (icd_code in icd_codes) {
      icd_presence_wide[, icd_code] <- icd_presence_wide[[icd_code]] | condition_presence_wide[[variable]]
    }
  }
}

# Filter to codes present in > 1,000 individuals

icd_cols_freq <- icd_presence_wide %>% colSums
icd_cols_1k_freq <- icd_cols_freq[icd_cols_freq >= 1000]

stop()

variables <- names(icd_cols_1k_freq) # colnames(icd_presence_wide)

icd_presence_wide <- icd_presence_wide[, variables]

condition_presence_table <- bind_cols(list(t2d_endpoint, icd_presence_wide))

# Add PRS
condition_presence_table <- add_t2d_scores_to_df(condition_presence_table, 'FINNGENID')$df

training_index <- !(k_fold_cv_indexes == fold_number)

result <- train_and_test_survival_prs_interaction(condition_presence_table, variables, training_index, alpha = 1)

saveRDS(result, paste0(save_path, "result_fold_", fold_number, ".rds"))
