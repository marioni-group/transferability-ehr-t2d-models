library(dplyr)

data_path <- "data/r11_condition/20231113/"

load_file <- function(filename) {
  if (endsWith(tolower(filename), ".rds")) {
    readRDS(paste0(data_path, filename))
  } else {
    stop(paste0("Error: unsupported file extension in load_file: ", filename))
  }
}

condition_preprocessing_result <- load_file("condition_preprocessing_result.rds")

condition_presence_wide <- condition_preprocessing_result$condition_presence_wide
t2d_endpoint <- condition_preprocessing_result$t2d_endpoint

t2d_endpoint$T2D_INCLAVO <- ifelse(t2d_endpoint$time_to_event > 5, 0, t2d_endpoint$T2D_INCLAVO)

condition_preprocessing_result <- NULL
gc()

incident_cases <- t2d_endpoint %>% filter(T2D_INCLAVO == 1)

cases_mean_tte <- mean(incident_cases$time_to_event)
cases_sd_tte <- sd(incident_cases$time_to_event)

cases_mean_onset_age <- mean(incident_cases$T2D_INCLAVO_AGE)
cases_sd_onset_age <- sd(incident_cases$T2D_INCLAVO_AGE)

cases_n_female <- sum(incident_cases$SEX == 1)
cases_percentage_female <- sum(incident_cases$SEX == 1)/nrow(incident_cases) * 100

censored <- t2d_endpoint %>% filter(T2D_INCLAVO == 0)

censored_mean_ttc <- mean(censored$time_to_event)
censored_sd_ttc <- sd(censored$time_to_event)

censored_mean_censoring_age <- mean(censored$T2D_INCLAVO_AGE)
censored_sd_censoring_age <- sd(censored$T2D_INCLAVO_AGE)

censored_n_female <- sum(censored$SEX == 1)
censored_percentage_female <- sum(censored$SEX == 1)/nrow(censored) * 100

df <- data.frame(
  Cases = c(
    n = nrow(incident_cases),
    mean_tte = cases_mean_tte,
    sd_tte = cases_sd_tte,
    mean_event_age = cases_mean_onset_age,
    sd_event_age = cases_sd_onset_age,
    n_female = cases_n_female,
    percentage_female = cases_percentage_female
  ),
  Censored = c(
    n = nrow(censored),
    mean_tte = censored_mean_ttc,
    sd_tte = censored_sd_ttc,
    mean_event_age = censored_mean_censoring_age,
    sd_event_age = censored_sd_censoring_age,
    n_female = censored_n_female,
    percentage_female = censored_percentage_female
  ))