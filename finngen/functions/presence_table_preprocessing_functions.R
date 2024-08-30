library(dplyr)
library(tidyr)
library(purrr)
library(glue)

long_to_wide <- function(df, person_id_table = NULL, frequency = FALSE, id_column = "concept_name") {
  if (is.null(person_id_table)) {
    person_id_table <- df %>% select(person_id) %>% distinct
  }
  
  # First group by concept. Then split into list of tbls
  group_split_result <- df %>% 
    group_by(across(all_of(id_column))) %>% 
    select(person_id, all_of(id_column)) 
  
  if (!frequency) {
    group_split_result <- group_split_result %>% distinct
  } 
  group_split_result <- group_split_result %>% group_split
  
  cids <- sapply(group_split_result, function(t) {
    name <- pull(t, all_of(id_column))[[1]]
    name
  })
  
  counts <- sapply(group_split_result, function(t) {
    nrow(t)
  })
  
  names(counts) <- cids

  result <- lapply(group_split_result, function(t) {
    name <- pull(t, all_of(id_column))[[1]]
    if (frequency) {
      t <- t %>% count(person_id) %>% rename('{{name}}':= n)
    } else {
      t <- t %>% mutate('{{name}}' := 1) %>% select(-all_of(id_column))
    }
    person_id_table %>% left_join(t, by = 'person_id') %>% select(-person_id)
  })
  
  result <- bind_cols(result)
  result[is.na(result)] <- 0
  result <- bind_cols(list(person_id_table, result))
  result <- rename_with(result, ~ gsub('"', '', .x)) # Remove quotes from column names
  list(table_wide = result, counts = counts)
}

presence_table_long_to_wide <- function(df, person_id_table = NULL, frequency = FALSE) {
  if (is.null(person_id_table)) {
    person_id_table <- df %>% select(person_id) %>% distinct
  }
  
  # First group by concept. Then split into list of tbls
  group_split_result <- df %>% 
    group_by(cid) %>% 
    select(person_id, cid) 
  if (!frequency) {
    group_split_result <- group_split_result %>% distinct
  } 
  group_split_result <- group_split_result %>% group_split

  cids <- sapply(group_split_result, function(t) {
    name <- pull(t, cid)[[1]]
    name
  })
  
  counts <- sapply(group_split_result, function(t) {
    nrow(t)
  })
  
  names(counts) <- cids
  
  result <- lapply(group_split_result, function(t) {
    name <- pull(t, cid)[[1]]
    if (frequency) {
      t <- t %>% count(person_id) %>% rename('{{name}}':= n)
    } else {
      t <- t %>% mutate('{{name}}' := 1) %>% select(-cid)
    }
    person_id_table %>% left_join(t, by = 'person_id') %>% select(-person_id)
  })

  result <- bind_cols(result)
  result <- bind_cols(list(person_id_table, result)) %>% replace(is.na(.), 0) 
  result <- rename_with(result, ~ gsub('"', '', .x)) # Remove quotes from column names
  list(table_wide = result, counts = counts)
}