library(dplyr)
library(Matrix)

df <- readRDS("../data/snomed_condition_frequency_table.rds")


# Requires topic model fold and cox lasso fold
generate_snomed_grouping_by_topic_table <- function(tm_fold_number, cl_fold_number) {
  snomed_to_topic_table <- read.csv(paste0("../models/topic_model_grouping_tm_all_folds/condition_topics_unique_fold_", tm_fold_number, ".csv"))

  snomed_to_topic_table$topic <- paste0("t", snomed_to_topic_table$topic)

  lasso_coefs <- read.csv(paste0("../models/topic_model_grouping_tm_all_folds/snomed_grouped_by_topic_lasso_coefs_result_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".csv"))

  df_topic <- data.frame(matrix(0, nrow=nrow(df), ncol = nrow(lasso_coefs)))
  colnames(df_topic) <- lasso_coefs %>% pull(X)

  snomed_to_topic <- function(x) {
    topic <- snomed_to_topic_table %>% filter(term == x) %>% pull(topic)
  }

  snomed_to_topic_lookup <- snomed_to_topic_table$topic
  names(snomed_to_topic_lookup) <- snomed_to_topic_table$term

  variables <- colnames(df)

  for (variable in variables) {
    topic <- ifelse(variable %in% names(snomed_to_topic_lookup), snomed_to_topic_lookup[[variable]], "")

    if (topic %in% colnames(df_topic)) {
      df_topic[, topic] <- df_topic[, topic] | (df[, variable] > 0)
    }
  }

  df_topic <- lapply(df_topic, as.numeric)

  # Add back in Cox table columns
  saveRDS(df_topic, paste0("../data/presence_snomed_grouped_by_topic_model_tm_all_folds_cl_fold_", cl_fold_number, "_tm_fold_", tm_fold_number, ".rds"))
}

for (tm_i in 1) {
  for (cl_i in 1:5) {
    generate_snomed_grouping_by_topic_table(tm_i, cl_i)
  }
}
