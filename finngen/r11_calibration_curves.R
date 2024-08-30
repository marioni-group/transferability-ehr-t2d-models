library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rms)

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

icd_results <- lapply(1:number_of_folds, function(i) {
  load_result(icd_grouping_result_path, i)
})
  
topic_model_results <- lapply(1:number_of_folds, function(i) {
  readRDS(paste0(topic_model_grouping_result_path, "result_cl_fold_", i, "_tm_fold_1.rds"))
})
  
null_calibration_curve <- CalibrationCurves::valProbggplot(snomed_results[[1]]$predictions$null_prediction, snomed_results[[1]]$test_surv$status,
                                                           statloc = FALSE,
                                                           CL.BT = TRUE,
                                                           CL.smooth = "fill",
                                                           d0lab = "",
                                                           d1lab = "",
                                                           legendloc = FALSE,
                                                           smooth = "loess")

ehr_calibration_curve <- CalibrationCurves::valProbggplot(snomed_results[[1]]$predictions$ehr_prediction, snomed_results[[1]]$test_surv$status,
                                                          statloc = FALSE,
                                                          CL.BT = TRUE,
                                                          CL.smooth = "fill",
                                                          d0lab = "",
                                                          d1lab = "",
                                                          legendloc = FALSE,
                                                          smooth = "loess")

icd_calibration_curve <- CalibrationCurves::valProbggplot(icd_results[[1]]$predictions$ehr_prediction, icd_results[[1]]$test_surv$status,
                                                          statloc = FALSE,
                                                          CL.BT = TRUE,
                                                          CL.smooth = "fill",
                                                          d0lab = "",
                                                          d1lab = "",
                                                          legendloc = FALSE,
                                                          smooth = "loess")

topic_calibration_curve <- CalibrationCurves::valProbggplot(topic_model_results[[1]]$predictions$ehr_prediction, topic_model_results[[1]]$test_surv$status,
                                                            statloc = FALSE,
                                                            CL.BT = TRUE,
                                                            CL.smooth = "fill",
                                                            d0lab = "",
                                                            d1lab = "",
                                                            legendloc = FALSE,
                                                            smooth = "loess")

calibration_plot <- ggarrange(null_calibration_curve$ggPlot + 
                                xlab(""), 
                              ehr_calibration_curve$ggPlot + 
                                xlab("") + 
                                ylab("") +
                                theme(legend.position = "none"), 
                              icd_calibration_curve$ggPlot +
                                theme(legend.position = "none"), 
                              topic_calibration_curve$ggPlot + 
                                ylab("") +
                                theme(legend.position = "none"))

pdf("plots/r11/publication_plots/calibration_plot.pdf")
print(calibration_plot)
dev.off()
