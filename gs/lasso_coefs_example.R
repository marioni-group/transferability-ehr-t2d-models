lasso_coefs <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/models/r11_model_files/1k_snomed/1k_snomed_cox_lasso_lasso_coefs_result_fold_1.csv")

gs_df <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table.rds")

lasso_coefs$present_in_gs <- sapply(lasso_coefs$X, function(x) {x %in% colnames(gs_df)})

write.csv(lasso_coefs, "/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/lasso_coefs_example.csv")

