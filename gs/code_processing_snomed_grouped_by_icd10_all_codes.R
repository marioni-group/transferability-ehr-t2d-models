library(dplyr)

df <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table.rds")

snomed_to_icd10_table <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/snomed_icd10_map.csv")
# Remove parts of term in brackets e.g. (finding)
snomed_to_icd10_table$term <- gsub("\\s*\\([^\\)]+\\)","",as.character(snomed_to_icd10_table$term))

generate_snomed_grouped_by_icd10_df <- function(fold_number) {
  lasso_coefs <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/icd_category_grouping/snomed_grouped_by_icd_1k_freq_cox_lasso_lasso_coefs_result_fold_", fold_number, ".csv"))

  df_icd10 <- data.frame(matrix(0, nrow=nrow(df), ncol = nrow(lasso_coefs)))
  colnames(df_icd10) <- lasso_coefs %>% pull(X)

  snomed_to_truncated_icd10 <- function(x) {
    icd10_string <- snomed_to_icd10_table %>% filter(term == x) %>% pull(truncated_icd10_code)
    icd10_string %>% strsplit(",") %>% unlist
  }

  snomed_to_icd10_lookup <- snomed_to_icd10_table$truncated_icd10_code
  names(snomed_to_icd10_lookup) <- snomed_to_icd10_table$term

  variables <- colnames(df)

  for (variable in variables) {
    icd_codes <- ifelse(variable %in% names(snomed_to_icd10_lookup), snomed_to_icd10_lookup[[variable]], "")
    icd_codes <- strsplit(icd_codes, split = ",")[[1]]

    if (length(icd_codes) > 0) {
      for (icd_code in icd_codes) {
        if (icd_code %in% colnames(df_icd10)) {
          df_icd10[, icd_code] <- df_icd10[[icd_code]] | df[[variable]]
        }
      }
    }
  }

  # Add back in Cox table columns
  df_icd10 <- cbind(df_icd10, df[, c("id", "Age", "sex", "Event", "time_to_event")])
  saveRDS(df_icd10, paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table_snomed_grouped_by_icd10_fold_", fold_number, ".rds"))
}
