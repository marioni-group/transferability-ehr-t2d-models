library(dplyr)
library(UpSetR)

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

fg_codes <- condition_presence_wide %>% colnames
# Remove "person_id"
fg_codes <- fg_codes[-1]

gs_codes <- readRDS("data/gs_unique_snomed_terms.rds")

upset_df <- fromList(list(GS = gs_codes, FG = fg_codes))
upset(upset_df, text.scale = 2)

fg_code_frequencies <- condition_presence_wide[, -1] %>% apply(2, sum)
fg_1k_code_frequencies <- fg_code_frequencies[fg_code_frequencies >= 1000]
fg_1k_codes <- fg_codes[fg_codes %in% names(fg_1k_code_frequencies)]

fg_1k_upset_df <- fromList(list(GS = gs_codes, FG_1k_freq = fg_1k_codes))
upset(fg_1k_upset_df, text.scale = 2)