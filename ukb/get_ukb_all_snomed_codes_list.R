condition_df <- readRDS("../data/snomed_condition_frequency_table.rds")
ukb_snomed_codes <- colnames(condition_df)[-1] # Remove "eid"
saveRDS(ukb_snomed_codes, "../results/ukb_all_snomed_codes_list.rds")
