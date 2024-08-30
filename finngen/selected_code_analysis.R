library(UpSetR)

freq_1k_counts <- lapply(1:5, function(i) {
  readRDS(paste0("data/r11_condition/20240723_selected_code_analysis/freq_1k_code_counts_fold_", i, ".rds"))
})



freq_1k_lasso_selected_codes <- lapply(1:5, function(i) {
  read.csv(paste0("models/r11_condition/20231113/1k_snomed_cox_lasso/extracted_model_tables/1k_snomed_cox_lasso_lasso_coefs_result_fold_", i, ".csv"))
})

code_to_freq <- function(codes, freqs) {
  sapply(codes, function(code) {
    freqs[[code]]
  })
}

freq_1k_code_freqs <- lapply(1:5, function(i) {
  code_to_freq(freq_1k_lasso_selected_codes[[i]]$X, freq_1k_counts[[i]])
})

freq_1k_codes <- lapply(freq_1k_counts, names)

names(freq_1k_codes) <- sapply(1:5, function(x) {paste0("fold_", x)})

lasso_selected_codes <- lapply(1:5, function(x) {
  freq_1k_lasso_selected_codes[[x]]$X
})
names(lasso_selected_codes) <- sapply(1:5, function(x) {paste0("fold_", x)})

pdf("plots/r11/publication_plots/fold_overlap_freq_1k.pdf", width = 5, height = 5)
print(upset(fromList(freq_1k_codes), order.by = "freq")) #, text.scale = 1))
dev.off()

pdf("plots/r11/publication_plots/lasso_fold_overlap.pdf", width = 5, height = 5)
print(upset(fromList(lasso_selected_codes), order.by = "freq")) # , text.scale = 1)
dev.off()

# Save code frequencies as csvs
lapply(1:5, function(fold_number) {
  df <- data.frame(freq_1k_counts[[fold_number]])
  write.csv(df, paste0("plots/r11/publication_plots/code_frequencies_fold_", fold_number, ".csv"))
})
