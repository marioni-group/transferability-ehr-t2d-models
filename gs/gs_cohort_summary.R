library(dplyr)

df <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/data/r11/gs_all_codes_presence_cox_table.rds")
df <- df %>% filter(time_to_event > 0)

df$AGE <- df$Age
df$SEX <- df$sex
t2d_prs <- read.csv("/Cluster_Filespace/Marioni_Group/Yipeng/INTERVENE/flagship/prs_output/T2D_PRS.sscore", sep = "\t")
t2d_prs$IID <- as.character(t2d_prs$IID)
# Add PRS
df <- df %>% left_join(t2d_prs, by = c("id" = "IID"))
df <- df %>% filter(!is.na(SCORE1_SUM))

cases <- df %>% filter(Event == 1)
controls <- df %>% filter(Event == 0)

tables <- list(cases = cases,
               controls = controls)

n <- sapply(tables, nrow)

mean_tte <- sapply(tables, function(table) {mean(table$time_to_event)})
sd_tte <- sapply(tables, function(table) {sd(table$time_to_event)})

mean_age <- sapply(tables, function(table) {mean(table$Age + table$time_to_event)})
sd_age <- sapply(tables, function(table) {sd(table$Age + table$time_to_event)})

n_female <- sapply(tables, function(table) {sum(table$sex == 1)})
percent_female <- sapply(tables, function(table) {sum(table$sex == 1) / nrow(table)})

summary_df <- data.frame(n = n,
                 mean_tte = mean_tte,
                 sd_tte = sd_tte,
                 mean_age = mean_age,
                 sd_age = sd_age,
                 n_female = n_female,
                 percent_female = percent_female)

write.csv(summary_df, "/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/ehr_t2d/results/r11/gs_cohort_summary.csv")
