library(dplyr)
library(lubridate)
library(zoo)

code_annotations <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/2023-11-09_codes_annotation.txt", sep = "\t")
primary <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/2023-11-09_primary_gp_annotated.txt", sep = "\t")
secondary <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/2023-11-09_secondary_smr01_annotated.txt", sep = "\t")

secondary <- secondary %>% rename(dt_ym = dt1_ym)

selected_cols <- c("id", "dt_ym", "concept_code")

combined <- rbind(primary[, selected_cols], secondary[, selected_cols])

snomed_concept_term_table <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/snomed_concept_term_table.csv")
# Remove parts of term in brackets e.g. (finding)
snomed_concept_term_table$term <- gsub("\\s*\\([^\\)]+\\)","",as.character(snomed_concept_term_table$term))
row.names(snomed_concept_term_table) <- as.character(snomed_concept_term_table$conceptId)

concept_term_lookup <- snomed_concept_term_table$term
names(concept_term_lookup) <- as.character(snomed_concept_term_table$conceptId)

# Replace concept_code with term
combined$concept_code <- as.character(combined$concept_code)
combined$concept_code <- sapply(combined$concept_code, function(x) {concept_term_lookup[[x]]})


unique_concept_ids <- combined %>% pull(concept_code) %>% unique

appt <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/2021-07-30_appt.csv")

appt$id <- as.character(appt$id)
combined$id <- as.character(combined$id)


# Filter codes to before appointment date
combined <- combined %>% left_join(appt, by = "id")
combined$dt_ym <- combined$dt_ym %>% as.character %>% as.yearmon(format = "%Y%m")
combined <- combined %>% filter(!is.na(dt_ym))
combined$pre_baseline <- combined$dt_ym < combined$appt
combined <- combined %>% filter(pre_baseline)

# Convert long condition table to wide presence table
person_id_table <- combined %>% select(id) %>% distinct

group_split_result <- combined %>%
  group_by(concept_code) %>%
  select(id, concept_code) %>% distinct %>% group_split

cids <- sapply(group_split_result, function(t) {
  name <- pull(t, concept_code)[[1]]
  name
})

counts <- sapply(group_split_result, function(t) {
  nrow(t)
})

names(counts) <- cids

result <- lapply(group_split_result, function(t) {
  name <- pull(t, concept_code)[[1]]

  t <- t %>% mutate('{{name}}' := 1) %>% select(-concept_code)
  person_id_table %>% left_join(t, by = "id") %>% select(-id)
})

result <- bind_cols(result)
result <- bind_cols(list(person_id_table, result)) %>% replace(is.na(.), 0)
result <- rename_with(result, ~ gsub('"', '', .x))

coxTable <- read.csv('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/using_methylpiper/src/20k/preprocessing/Diabetes_all_cases_150322.csv')
coxTable$Sample_Name <- as.character(coxTable$Sample_Name)
result <- result %>% left_join(coxTable[,c("Sample_Name", "Age", "sex", "Event", "time_to_event", "dead")], by = c("id" = "Sample_Name"))
result <- result %>% filter(!is.na(Event))
result$sex <- as.numeric(result$sex == "F")
result$Event_before_thresholding <- result$Event
result$Event <- ifelse(result$time_to_event > 5, 0, result$Event)

saveRDS(result, "/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table_with_dead_and_event_before_thresholding_columns.rds")
