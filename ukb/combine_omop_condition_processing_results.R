library(dplyr)
library(Matrix)

eid_ranges <- c(
    "1-50000",
    "50001-1e+05",
    "100001-150000",
    "150001-2e+05",
    "200001-250000",
    "250001-3e+05",
    "300001-350000",
    "350001-4e+05",
    "400001-450000",
    "450001-5e+05",
    "500001-502359"
)

results <- lapply(eid_ranges, function(eid_range) {
    readRDS(paste0("../data/condition_frequency_table_eid_index_", eid_range, "_with_dummy.rds"))[-1, ] # Remove dummy row
})


condition_frequency_df <- do.call(rbind, results)

covariates_df <- readRDS("../data/covariates.rds")
cox_df <- readRDS("../data/t2d_cox_table.rds")


row.names(condition_frequency_df) <- condition_frequency_df[, "eid"]
condition_frequency_df <- condition_frequency_df[as.character(cox_df$eid), ]

saveRDS(condition_frequency_df, "../data/snomed_condition_frequency_table.rds")
