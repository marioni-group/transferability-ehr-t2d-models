library(dplyr)
library(ggplot2)

snomed_result_path <- "models/r11_condition/20231113/1k_snomed_cox_lasso/"
icd_grouping_result_path <- "models/r11_condition/20231115/snomed_grouped_by_icd_1k_freq_cox_lasso/"
topic_model_grouping_result_path <- "models/r11_condition/20231124/cox_lasso_models/"

number_of_folds <- 5

load_result <- function(folder_path, fold_number) {
  readRDS(paste0(folder_path, "result_fold_", fold_number, ".rds"))
}

snomed_results <- lapply(1:number_of_folds, function(i) {
  load_result(snomed_result_path, i)
})


topic_model_grouping_results <- lapply(1:3, function(tm_fold_number) {
  lapply(1:number_of_folds, function(i) {
    readRDS(paste0(topic_model_grouping_result_path, "result_cl_fold_", i, "_tm_fold_", tm_fold_number, ".rds"))
  })
})

icd_grouping_results <- lapply(1:number_of_folds, function(i) {
  load_result(icd_grouping_result_path, i)
})

null_aucs <- sapply(snomed_results, function(x) {
  MLmetrics::PRAUC(x$rocs$null$predictor, x$rocs$null$response)
})

ehr_snomed_1k_aucs <- sapply(snomed_results, function(x) {
  MLmetrics::PRAUC(x$rocs$ehr$predictor, x$rocs$ehr$response)
})
ehr_snomed_1k_auc_df <- data.frame(Model = "SNOMED 1k Frequency", PRAUC = ehr_snomed_1k_aucs)

ehr_icd_grouping <- sapply(icd_grouping_results, function(x) {
  MLmetrics::PRAUC(x$rocs$ehr$predictor, x$rocs$ehr$response)
})
ehr_icd_grouping_auc_df <- data.frame(Model = "ICD Category Grouping", PRAUC = ehr_icd_grouping)

ehr_topic_grouping <- sapply(topic_model_grouping_results, function(result_by_tm_fold) {
  sapply(result_by_tm_fold, function(x) {
    MLmetrics::PRAUC(x$rocs$ehr$predictor, x$rocs$ehr$response)
  })
}) %>% as.vector
ehr_topic_grouping_auc_df <- data.frame(Model = "Topic Grouping", PRAUC = ehr_topic_grouping)

auc_df <- rbind(ehr_snomed_1k_auc_df,
                ehr_icd_grouping_auc_df,
                ehr_topic_grouping_auc_df)

auc_df$Model <- factor(auc_df$Model, levels = c("SNOMED 1k Frequency",
                                                "ICD Category Grouping",
                                                "Topic Grouping"))

auc_plot <- auc_df %>% ggplot(aes(x = Model, y = PRAUC)) +
  geom_point(alpha = 0.5) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 22,
    fill = "red"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_minimal(base_size = 20) +
  xlab("Model (all include Age + Sex as covariates)") +
  ylab("PRAUC") +
  ggtitle("PRAUC in FinnGen")

pdf("plots/r11/publication_plots/finngen_prauc_plot.pdf")
print(auc_plot)
dev.off()

auc_means <- auc_df %>% group_by(Model) %>% summarise(Mean = mean(PRAUC))
prs_aucs <- sapply(snomed_results, function(x) {
  MLmetrics::PRAUC(x$rocs$prs$predictor, x$rocs$prs$response)
})
prs_aucs_mean <- mean(prs_aucs)