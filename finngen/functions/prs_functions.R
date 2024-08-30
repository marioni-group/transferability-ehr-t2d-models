library(readr)
library(dplyr)
library(tidyr)

add_t2d_scores_to_df <- function(df, person_id_col, finngen_data_freeze = "R11") {
  scores <- read_tsv(paste0('/finngen/library-red/finngen_', finngen_data_freeze, '/prs_1.0/data/finngen_', finngen_data_freeze, '_Mahajan.NatGenet2018b.T2D.European.txt.sscore'))
  # setNames used here so person_id_col can be used as a variable to specify join column
  df <- df %>% left_join(scores %>% select(IID, SCORE1_AVG), by = setNames('IID', person_id_col))
  list(scores_table = scores, df = df)
}

