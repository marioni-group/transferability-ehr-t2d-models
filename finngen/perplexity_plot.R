library(ggplot2)
library(dplyr)

perplexity_plot_dfs <- lapply(1:3, function(i) {
  p <- readRDS(paste0("models/r11_condition/20231124/topic_models/perplexity_plot_df_fold_", i, ".rds"))
  p$cv_iter <- i
  p
})

perplexity_plot_df <- do.call(rbind, perplexity_plot_dfs)

perplexity_plot <- perplexity_plot_df %>% ggplot(aes(x = n_topics, y = Perplexity, color = as.factor(cv_iter))) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal(base_size = 20) +
  xlab("Number of topics") + 
  ylab("Perplexity") +
  guides(color = guide_legend(title = "CV Iteration"))

pdf("plots/r11/publication_plots/perplexity_plot.pdf")
print(perplexity_plot)
dev.off()