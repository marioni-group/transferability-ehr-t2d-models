library(dplyr)
library(tidyr)
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

topic_model_grouping_results <- lapply(1:number_of_folds, function(i) {
  readRDS(paste0(topic_model_grouping_result_path, "result_cl_fold_", i, "_tm_fold_1.rds"))
})

icd_grouping_results <- lapply(1:number_of_folds, function(i) {
  load_result(icd_grouping_result_path, i)
})

snomed_result_aucs <- sapply(snomed_results, function(x) {x$aucs}) %>% t %>% as.data.frame
icd_grouping_result_aucs <- sapply(icd_grouping_results, function(x) {x$aucs}) %>% t %>% as.data.frame
topic_model_grouping_result_aucs <- sapply(topic_model_grouping_results, function(x) {x$aucs}) %>% t %>% as.data.frame

auc_df <- data.frame(`Null na` = snomed_result_aucs$null_auc %>% unlist,
                     `PRS na` = snomed_result_aucs$prs_auc %>% unlist,
                     `EHR Original` = snomed_result_aucs$ehr_auc %>% unlist,
                     `EHR ICD` = icd_grouping_result_aucs$ehr_auc %>% unlist,
                     `EHR Topic` = topic_model_grouping_result_aucs$ehr_auc %>% unlist,
                     `Full Original` = snomed_result_aucs$full_auc %>% unlist,
                     `Full ICD` = icd_grouping_result_aucs$full_auc %>% unlist,
                     `Full Topic` = topic_model_grouping_result_aucs$full_auc %>% unlist,
                     `Interaction Original` = snomed_result_aucs$interaction_auc %>% unlist,
                     `Interaction ICD` = icd_grouping_result_aucs$interaction_auc %>% unlist,
                     `Interaction Topic` = topic_model_grouping_result_aucs$interaction_auc %>% unlist)

auc_df$fold <- 1:5
auc_plot_df <- auc_df %>% gather(key = "variable", value = "value", -c(fold)) %>% separate(variable, c("type", "grouping"))

type_to_model_string <- function(x) {
  sapply(x, function(type) {
    if (type == "Null") {
      "Age + Sex only"
    } else if (type == "Full") {
      "PRS + EHR"
    } else if (type == "Interaction") {
      "PRS * EHR"
    } else if (type == "EHR") {
      "EHR"
    } else if (type == "PRS") {
      "PRS"
    }
  })
}

auc_plot_df <- auc_plot_df %>% mutate(type = type_to_model_string(type)) %>% mutate(grouping = factor(grouping, levels = c("na", "Original", "ICD", "Topic")))

auc_plot_df

auc_plot_df <- auc_plot_df %>% 
  mutate(Model = paste0(auc_plot_df$type, "_", auc_plot_df$grouping))

auc_plot_df$Model <- sapply(auc_plot_df$Model, function(model_string) {
  if (model_string == "Age + Sex only_na") {
    "Age + Sex only"
  } else if (model_string == "PRS_na") {
    "PRS"
  } else {
    model_string
  }
})
auc_plot_df

auc_plot <- auc_plot_df %>% ggplot(aes(x = as.factor(Model), y = value)) + 
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
  ylab("AUC") +
  ggtitle("AUCs in FinnGen")
  

auc_plot_df_ehr <- auc_plot_df %>% 
  filter(Model %in% c("Age + Sex only", "EHR_Original", "EHR_ICD", "EHR_Topic")) %>% 
  mutate(Model = sapply(Model, function(x) {
    if (x == "Age + Sex only") {
      "Age + Sex only"
    } else if (x == "EHR_Original") {
      "SNOMED 1k Frequency"
    } else if (x == "EHR_ICD") {
      "ICD Category Grouping"
    } else if (x == "EHR_Topic") {
      "Topic Grouping"
    }
  }))

auc_plot_ehr <- auc_plot_df_ehr %>% ggplot(aes(x = as.factor(Model), y = value)) + 
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
  ylab("AUC") +
  ggtitle("AUCs in FinnGen")

pdf("plots/r11/publication_plots/finngen_auc_plot.pdf")
print(auc_plot_ehr)
dev.off()