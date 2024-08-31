library(dplyr)

covariates_df <- readRDS("../data/covariates.rds")
t2d_cox_df <- readRDS("../data/t2d_cox_table.rds")

df <- t2d_cox_df %>% left_join(covariates_df)

# Threshold cases to tte <= 5
df$Event <- ifelse(df$time_to_event <= 5, df$Event, 0)

df$event_age <- df$age + df$time_to_event

cases <- df %>% filter(Event == 1)
controls <- df %>% filter(Event == 0)

calculate_stats <- function(x) {
    list(
        n = nrow(x),
        mean_tte = mean(x$time_to_event),
        sd_tte = sd(x$time_to_event),
        mean_age = mean(x$age),
        sd_age = sd(x$age),
        mean_event_age = mean(x$event_age),
        sd_event_age = sd(x$event_age),
        n_female = sum(x$sex == 1),
        percent_female = sum(x$sex == 1)/nrow(x)
    )
}

case_stats <- calculate_stats(cases)
control_stats <- calculate_stats(controls)

cohort_details_df <- rbind(data.frame(case_stats), data.frame(control_stats))
row.names(cohort_details_df) <- c("cases", "controls")

write.csv(cohort_details_df, "../results/ukb_cohort_details.csv")
