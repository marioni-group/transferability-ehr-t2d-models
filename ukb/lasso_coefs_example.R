lasso_coefs <- read.csv("../models/1k_snomed/1k_snomed_cox_lasso_lasso_coefs_result_fold_1.csv")

ukb_df <- readRDS("../data/snomed_condition_frequency_table.rds")

lasso_coefs$present_in_ukb <- sapply(lasso_coefs$X, function(x) {x %in% colnames(ukb_df)})

write.csv(lasso_coefs, "../data/lasso_coefs_example.csv")
