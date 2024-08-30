library(dplyr)

condition_codes <- readRDS("data/r11_condition/20231113/condition_preprocessing_result.rds")$condition_id_name_table
gc()
snomed_to_icd10_table <- read.csv("data/snomed_icd10_map.csv")

# Remove (disorder) and (finding) suffixes from term column
snomed_to_icd10_table[, "term"] <- sapply(snomed_to_icd10_table[, "term"], function(x) {
  x <- gsub(" \\(disorder\\)", "", x)
  gsub(" \\(finding\\)", "", x)
})

snomed_codes <- condition_codes %>% pull(concept_name)


snomed_to_truncated_icd10 <- function(x) {
  icd10_string <- snomed_to_icd10_table %>% filter(conceptId == x) %>% pull(truncated_icd10_code)
  icd10_string %>% strsplit(",") %>% unlist
}

snomed_to_icd10_lookup <- snomed_to_icd10_table$truncated_icd10_code
names(snomed_to_icd10_lookup) <- snomed_to_icd10_table$term

snomed_to_icd10 <- function(x) {
  lookup_result <- ifelse(x %in% names(snomed_to_icd10_lookup), snomed_to_icd10_lookup[[x]], "")
  # separate lookup results that have multiple ICD10 codes into a vector
  lookup_result %>% strsplit(",") %>% unlist
}

mapped_icd10 <- sapply(snomed_codes, snomed_to_icd10)

mapped_icd10_counts <- unlist(mapped_icd10) %>% table
unique_mapped_icd10 <- unique(unlist(mapped_icd10))

