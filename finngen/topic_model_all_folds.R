library(topicmodels)
library(dplyr)
library(readr)
library(ggplot2)
library(tidytext)
source("src/functions/util_functions.R")

n_topics <- 200

set.seed(42)
condition_preprocessing_result <- readRDS("data/r11_condition/20231113/condition_preprocessing_result.rds")
condition_frequency_wide <- condition_preprocessing_result$condition_frequency_wide
t2d_endpoint <- readRDS("data/r11_condition/20231113/t2d_endpoint_with_topic_model_cv_indexes.rds")
condition_preprocessing_result <- NULL
gc()

# Filter condition frequency table to topic modelling set and training folds
model_fitting_set <- condition_frequency_wide[t2d_endpoint$prediction_model_set_index == 0, ]

t2d_endpoint <- t2d_endpoint %>% filter(prediction_model_set_index == 0)

condition_frequency_wide <- NULL
gc()

condition_presence_counts <- ((model_fitting_set %>% select(-person_id)) > 0) %>% apply(2, sum)
condition_presence_1k <- condition_presence_counts[condition_presence_counts >= 100] %>% names

person_id_table_mfs <- model_fitting_set %>% select(person_id)

model_fitting_set <- model_fitting_set[, colnames(model_fitting_set) %in% condition_presence_1k]

gc()

non_zero_index_mfs <- model_fitting_set %>% apply(1, function(x) {any(x > 0)})

model_fitting_set <- model_fitting_set[non_zero_index_mfs, ]
gc()


print(paste0("Number of topics: ", n_topics))
print(system.time({tm <- LDA(model_fitting_set, k = n_topics, control = list(seed = 42))}))

tm <- readRDS("models/r11_condition/20240814_topic_model_all_folds/topic_model_200_topics_all_folds.rds")

process_topic_model <- function(topic_model, n_top_terms = 5) {
  
  condition_topics <- tidy(topic_model, matrix = "beta")
  
  condition_topics_max_beta <- condition_topics %>% group_by(term) %>% filter(beta == max(beta)) %>% add_count(name = "term_topic_count") %>% ungroup()
  
  # Some terms have the same probability for all topics and therefore all their max values are the same
  # Filter condition_topics_max_beta to only terms that have a unique max
  condition_topics_unique <- condition_topics_max_beta %>% filter(term_topic_count == 1)
  
  ap_top_terms <- condition_topics %>%
    group_by(topic) %>%
    slice_max(beta, n = n_top_terms) %>% 
    ungroup() %>%
    arrange(topic, -beta)
  list(condition_topics = condition_topics, 
       condition_topics_max_beta = condition_topics_max_beta,
       condition_topics_unique = condition_topics_unique,
       ap_top_terms = ap_top_terms)
}

process_topic_model_result <- process_topic_model(tm)

saveRDS(process_topic_model_result, paste0("models/r11_condition/20240814_topic_model_all_folds/process_topic_model_result_", n_topics, "_topics_all_folds.rds"))
