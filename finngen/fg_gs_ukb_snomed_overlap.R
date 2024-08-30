library(gplots)
library(dplyr)

data_path <- "data/r11_condition/20231113/"

load_file <- function(filename) {
  if (endsWith(tolower(filename), ".rds")) {
    readRDS(paste0(data_path, filename))
  } else {
    stop(paste0("Error: unsupported file extension in load_file: ", filename))
  }
}

condition_preprocessing_result <- load_file("condition_preprocessing_result.rds")

condition_presence_wide <- condition_preprocessing_result$condition_presence_wide

condition_preprocessing_result <- NULL
gc()

fg_snomed_codes <- condition_presence_wide %>% 
  select(-person_id) %>% 
  colnames

gs_snomed_codes <- readRDS("data/r11_condition/gs_all_snomed_codes_list.rds")
ukb_snomed_codes <- readRDS("data/r11_condition/ukb_all_snomed_codes_list.rds")

venn_table <- list(FinnGen = fg_snomed_codes,
                   GS = gs_snomed_codes,
                   UKB = ukb_snomed_codes) %>%
  venn

