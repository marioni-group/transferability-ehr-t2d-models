library(dplyr)
library(readr)
library(bigrquery)

create_condition_presence_table <- function(baseline, lookback=NULL, frequency = FALSE, snomed_only = TRUE) {
  bq_auth(scopes='https://www.googleapis.com/auth/bigquery')
  
  projectid = 'fg-production-sandbox-6'
  
  # If frequency is set to TRUE, don't include DISTINCT in the SQL SELECT
  distinct_string <- ifelse(frequency, "", "DISTINCT")
  
  if (is.null(lookback)) {
    sql <- paste0("SELECT ", distinct_string, " person_id, condition_concept_id, concept_name, concept_code, vocabulary_id FROM
            (SELECT person_id, condition_concept_id, concept_name, concept_code, condition_start_DATETIME, vocabulary_id
            FROM finngen-production-library.finngen_omop_r11.condition_occurrence
            LEFT JOIN finngen-production-library.finngen_omop_r11.concept
            ON condition_concept_id = concept_id
            WHERE condition_start_DATETIME < '", baseline, "')
            LIMIT 50000000")
  } else {
    sql <- paste0("SELECT ", distinct_string, " person_id, condition_concept_id, concept_name, concept_code, vocabulary_id FROM
            (SELECT person_id, condition_concept_id, concept_name, concept_code, condition_start_DATETIME, vocabulary_id
            FROM finngen-production-library.finngen_omop_r11.condition_occurrence
            LEFT JOIN finngen-production-library.finngen_omop_r11.concept
            ON condition_concept_id = concept_id
            WHERE condition_start_DATETIME < '", baseline, "' AND condition_start_DATETIME >= '", lookback, "')
            LIMIT 50000000")
  }
  tb <- bq_project_query(projectid, sql)
  df <- bq_table_download(tb, page_size = 5000)
  
  if (snomed_only) {
    df <- df %>% filter(vocabulary_id == "SNOMED")
  }
  
  list(tb = tb, df = df)
}
