library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

library(bigrquery)
library(UpSetR)

source('src/functions/r11_cot_preprocessing_functions.R')

Sys.setenv(VROOM_CONNECTION_SIZE = 131072*8)


set.seed(42)

# Redesigned function using BigQuery rather than reading in whole endpoint table
get_endpoint_with_person <- function(endpoint_name, person_table = NULL) {
  if (is.null(person_table)) {
    person_table <- read_tsv("/finngen/library-red/finngen_R11/phenotype_1.0/data/finngen_omop/finngen_R11_person.csv")
  }
  
  
  projectid <- "fg-production-sandbox-6"
  options(gargle_oauth_cache = FALSE)
  bq_auth(scopes = "https://www.googleapis.com/auth/bigquery")
  

  sql<-paste0("
    SELECT *
    FROM finngen-production-library.sandbox_tools_r11.endpoint_cohorts_r11_v1
    WHERE ENDPOINT = '", endpoint_name, "'
    LIMIT 10000000
  ")
  
  tb <- bq_project_query(projectid, sql)
  df <- bq_table_download(tb)
  
  
  
  additional_columns_sql <- "
    SELECT FINNGENID, BL_YEAR, BL_AGE, SEX, HEIGHT, HEIGHT_AGE, WEIGHT, WEIGHT_AGE, FU_END_AGE, DEATH, DEATH_AGE, DEATH_YEAR, APPROX_BIRTH_DATE
    FROM finngen-production-library.sandbox_tools_r11.finngenid_info_r11_v1
    LIMIT 10000000
  "
  
  additional_columns_tb <- bq_project_query(projectid, additional_columns_sql)
  additional_columns_df <- bq_table_download(additional_columns_tb)
  
  new_age_colname <- paste0(endpoint_name, "_AGE")
  new_year_colname <- paste0(endpoint_name, "_YEAR")
  new_nevt_colname <- paste0(endpoint_name, "_NEVT")
  
  df <- df %>% rename_with(~endpoint_name, CONTROL_CASE_EXCL) %>%
    rename_with(~new_age_colname, AGE) %>%
    rename_with(~new_year_colname, YEAR) %>%
    rename_with(~new_nevt_colname, NEVT)
  
  df <- df %>% inner_join(additional_columns_df)
  
  # Reorder person_table so person_source_value matches FINNGENID in df
  person_table_index <- match(df$FINNGENID, person_table$person_source_value)
  person_table <- person_table[person_table_index,]
  
  df <- cbind(df, person_table)
  
  list(endpoint_with_person = df,
       person_table = person_table)
}

# Used in get_endpoint_filtered_by_baseline
calculate_onset_or_censoring_datetime <- function(df, age_colname = "T2D_INCLAVO_AGE") {
  df <- df %>% mutate(onset_or_censoring_datetime = birth_datetime + duration(.data[[age_colname]], units = 'years'))
  df
}

# Used in get_endpoint_filtered_by_baseline
filter_onset_or_censoring_datetime_before_baseline <- function(df, baseline = '2014-02-01', end_date = NULL) {
  df <- df %>% filter(onset_or_censoring_datetime > baseline)
  if (!is.null(end_date)) {
    df <- df %>% filter(onset_or_censoring_datetime < end_date)
  }
  df
}

# Returns the input endpoint table with individuals filtered out if T2D onset or censoring datetime < baseline
get_endpoint_filtered_by_baseline <- function(endpoint_with_person, baseline = '2014-02-01', end_date = NULL) {
  endpoint_with_onset_or_censoring_datetime <- calculate_onset_or_censoring_datetime(endpoint_with_person)
  endpoint_onset_or_censoring_filtered <- filter_onset_or_censoring_datetime_before_baseline(endpoint_with_onset_or_censoring_datetime, baseline = baseline, end_date = end_date)
  endpoint_onset_or_censoring_filtered
}

filter_to_alive <- function(df, threshold_date) {
  death_table <- read_tsv('/finngen/library-red/finngen_R11/phenotype_1.0/data/finngen_omop/finngen_R11_death.csv')

  df <- df %>% left_join(death_table %>% select(person_id, death_date), by = 'person_id') %>% filter(is.na(death_date) | (death_date >= threshold_date)) %>% select(-death_date)
  
  df
}

# For removing people with any of the specified endpoints before baseline
remove_additional_endpoints_before_baseline <- function(endpoint_with_person, additional_endpoint_names, person_table = NULL, baseline = '2014-02-01', end_date = NULL) {
  additional_endpoint_age_colname <- additional_endpoint_names %>% sapply(function(x) {paste0(x, '_AGE')})
  names(additional_endpoint_age_colname) <- NULL
  additional_endpoint_year_colname <- additional_endpoint_names %>% sapply(function(x) {paste0(x, '_YEAR')})
  names(additional_endpoint_year_colname) <- NULL
  additional_endpoint_nevt_colname <- additional_endpoint_names %>% sapply(function(x) {paste0(x, '_NEVT')})
  names(additional_endpoint_nevt_colname) <- NULL
  
  additional_endpoint_tables <- lapply(1:length(additional_endpoint_names), function(i) {
    current_additional_endpoint <- additional_endpoint_names[[i]]
    current_age_colname <- additional_endpoint_age_colname[[i]]
    current_year_colname <- additional_endpoint_year_colname[[i]]
    current_nevt_colname <- additional_endpoint_nevt_colname[[i]]
    
    endpoint_table <- get_endpoint_with_person(current_additional_endpoint)$endpoint_with_person %>%
      select(all_of(c("FINNGENID", current_additional_endpoint, current_age_colname, current_year_colname, current_nevt_colname)))
  })
  
  additional_endpoint_table <- Reduce(inner_join, additional_endpoint_tables)
  
  if (is.null(person_table)) {
    person_table <- read_tsv('/finngen/library-red/finngen_R11/phenotype_1.0/data/finngen_omop/finngen_R11_person.csv')
  }
  
  additional_endpoint_table <- additional_endpoint_table %>% left_join(person_table, by = c('FINNGENID' = 'person_source_value'))
  
  cases <- additional_endpoint_table %>% filter(rowSums(additional_endpoint_table %>% select(all_of(additional_endpoint_names)), na.rm = TRUE) > 0)
  

  # Calculate event datetime
  for (aen in additional_endpoint_names) {
    cases[[paste0(aen, '_datetime')]] <- cases$birth_datetime + duration(cases[[paste0(aen, '_AGE')]], units = 'years')
  }

  diabetes_datetime_columns <- sapply(additional_endpoint_names, function(x) {paste0(x, '_datetime')})
  names(diabetes_datetime_columns) <- NULL

  first_diabetes_datetime <- cases %>% select(all_of(diabetes_datetime_columns)) %>% apply(1, min)

  cases <- cases %>% mutate(first_diabetes_datetime = first_diabetes_datetime)

  cases_before_baseline <- cases %>% filter(first_diabetes_datetime < baseline)

  person_ids_to_remove <- cases_before_baseline %>% pull(person_id)

  endpoint_additional_removed <- endpoint_with_person %>% filter(!(person_id %in% person_ids_to_remove))

  list(cases = cases, cases_before_baseline = cases_before_baseline, first_diabetes_datetime = first_diabetes_datetime, endpoint_additional_removed = endpoint_additional_removed)
  }