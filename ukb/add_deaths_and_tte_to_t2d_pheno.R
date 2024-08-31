library(dplyr)
library(data.table)
library(lubridate)

phenotypes_data_path <- "../../../data/base_data/ukbb_78537/phenotypes/"

translations <- read.csv(paste0(phenotypes_data_path, "ukbb_78537_translations.tsv"), sep = "\t")

base_phenotypes_head <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), nrows = 5)

base_phenotypes_colnames <- colnames(base_phenotypes_head)

# Filter translations to those corresponding to columns that are present in base_phenotypes
translations_present <- translations %>% filter(UDI %in% base_phenotypes_colnames)

death_date_udis <- translations_present %>% filter(Description == "Date of death") %>% pull(UDI)

death_date_table <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = death_date_udis) %>% as.data.frame

death_date_table <- lapply(death_date_table, as.Date) %>% as.data.frame

eids <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = "eid")

death_date_table$eid <- eids$eid

baseline_date <- as.Date("2014-02-01")

death_date_table <- death_date_table[!is.na(death_date_table[, 1]), ]

death_before_baseline <- (death_date_table[, 1] <= baseline_date)

date_and_code_df <- readRDS("../data/t2d_date_and_code_df.rds")

# Remove participants with death before baseline
death_before_baseline_eids <- death_date_table[death_before_baseline, ]$eid

date_and_code_df <- date_and_code_df %>% filter(!(eid %in% death_before_baseline_eids))

# Add Event column
date_and_code_df$Event <- as.numeric(!is.na(date_and_code_df$date))

# Add time_to_event
death_date_table$death_date <- death_date_table[, 1]
date_and_code_df <- date_and_code_df %>% left_join(death_date_table %>% select(eid, death_date))

date_and_code_df$event_date <- sapply(1:nrow(date_and_code_df), function(i) {
    row <- date_and_code_df[i, ]
    if (row$Event) {
        row$date
    } else if (!is.na(row$death_date)) {
        row$death_date
    } else {
        as.Date("2022-10-31")
    }
}) %>% as.Date

date_and_code_df$time_to_event <- time_length(interval(as.Date("2014-02-01"), date_and_code_df$event_date), "years")

date_and_code_df <- date_and_code_df %>% filter(time_to_event > 0)

saveRDS(date_and_code_df, "../data/t2d_cox_table.rds")
