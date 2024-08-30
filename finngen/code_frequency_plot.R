library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(glmnet)
library(randomForest)
source('src/functions/prs_functions.R')
source('src/functions/model_fitting_functions.R')
source('src/functions/util_functions.R')

data_path <- "data/r11_condition/20231113/"
save_path <- "models/r11_condition/20240611/1k_snomed_cox_lasso/"
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

condition_counts <- condition_presence_wide %>% select(-person_id) %>% apply(2, sum) %>% sort(decreasing = TRUE)
condition_count_df <- data.frame(x = 1:length(condition_counts), frequency = condition_counts)

freq_plot <- condition_count_df %>%
  filter(frequency >= 5) %>%
  ggplot(aes(x = x, y = frequency)) +
    geom_point(alpha = 0.25, size = 2) +
    scale_x_continuous(breaks = seq(0, length(condition_counts), 1000)) +
    scale_y_continuous(breaks = seq(0, max(condition_counts), 2000)) + 
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + # , axis.ticks.x = element_blank()) +
    xlab("Codes sorted by frequency") +
    ylab("Frequency")
freq_plot


condition_count_df$log_frequency <- log(condition_count_df$frequency)

pdf("plots/r11/publication_plots/code_frequency_plot.pdf", width = 5, height = 6)
print(freq_plot)
dev.off()