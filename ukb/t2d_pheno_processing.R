library(dplyr)
library(data.table)

phenotypes_data_path <- "../../../data/base_data/ukbb_78537/phenotypes/"

translations <- read.csv(paste0(phenotypes_data_path, "ukbb_78537_translations.tsv"), sep = "\t")

base_phenotypes_head <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), nrows = 5)

base_phenotypes_colnames <- colnames(base_phenotypes_head)

# Filter translations to those corresponding to columns that are present in base_phenotypes
translations_present <- translations %>% filter(UDI %in% base_phenotypes_colnames)

main_diagnoses_descriptions <- c(
    "Diagnoses - ICD10Uses data-coding 19",
    "Diagnoses - ICD9Uses data-coding 87",
    "Date of first in-patient diagnosis - ICD10",
    "Date of first in-patient diagnosis - ICD9"
)

main_diagnoses_udis <- lapply(main_diagnoses_descriptions, function(x) {
    translations_present %>% filter(Description %in% x) %>% pull(UDI)
})

names(main_diagnoses_udis) <- c(
    "icd_10_code",
    "icd_9_code",
    "icd_10_date",
    "icd_9_date"
)

eids <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = "eid")

main_diagnoses_phenos <- lapply(main_diagnoses_udis, function(x) {
    gc()
    fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = x) %>% as.matrix
})

icd_10_e11_present <- apply(main_diagnoses_phenos$icd_10_code, c(1, 2), function(x) {
    startsWith(x, "E11")
})

e11_codes <- main_diagnoses_phenos$icd_10_code
e11_codes[!icd_10_e11_present] <- NA

e11_dates <- main_diagnoses_phenos$icd_10_date
e11_dates[!icd_10_e11_present] <- NA
e11_dates <- as.data.frame(e11_dates)
gc()
e11_dates <- lapply(e11_dates, as.Date) %>% as.data.frame
gc()

# For each row, find which column the min date is in
e11_min_date_indexes <- apply(e11_dates, 1, function(x) {
    ifelse(
        all(is.na(x)),
        NA,
        which(x == min(x, na.rm = TRUE))
    )
})

e11_date_first <- sapply(1:length(e11_min_date_indexes), function(i) {
    if (is.na(e11_min_date_indexes[[i]])) {
        NA
    } else {
        e11_dates[i, e11_min_date_indexes[[i]]]
    }
}) %>% as.Date

e11_codes_first <- sapply(1:length(e11_min_date_indexes), function(i) {
    if (is.na(e11_min_date_indexes[[i]])) {
        NA
    } else {
        e11_codes[i, e11_min_date_indexes[[i]]]
    }
})

icd_9_250_present <- apply(main_diagnoses_phenos$icd_9_code, c(1, 2), function(x) {
    # startsWith(x, "250")
    grepl("250[0-9]0", x) | grepl("250[0-9]2", x)
})

icd_9_250_codes <- main_diagnoses_phenos$icd_9_code
icd_9_250_codes[!icd_9_250_present] <- NA

icd_9_250_dates <- main_diagnoses_phenos$icd_9_date
icd_9_250_dates[!icd_9_250_present] <- NA
icd_9_250_dates <- as.data.frame(icd_9_250_dates)
gc()
icd_9_250_dates <- lapply(icd_9_250_dates, as.Date) %>% as.data.frame
gc()

# For each row, find which column the min date is in
icd_9_250_min_date_indexes <- apply(icd_9_250_dates, 1, function(x) {
    ifelse(
        all(is.na(x)),
        NA,
        which(x == min(x, na.rm = TRUE))
    )
})

icd_9_250_date_first <- sapply(1:length(icd_9_250_min_date_indexes), function(i) {
    if (is.na(icd_9_250_min_date_indexes[[i]])) {
        NA
    } else {
        icd_9_250_dates[i, icd_9_250_min_date_indexes[[i]]]
    }
}) %>% as.Date

icd_9_250_codes_first <- sapply(1:length(icd_9_250_min_date_indexes), function(i) {
    if (is.na(icd_9_250_min_date_indexes[[i]])) {
        NA
    } else {
        icd_9_250_codes[i, icd_9_250_min_date_indexes[[i]]]
    }
})


# Combine ICD-9 and ICD-10 to find min out of the two
icd_dates <- data.frame(icd9 = icd_9_250_date_first, icd10 = e11_date_first)
icd_codes <- data.frame(icd9 = icd_9_250_codes_first, icd10 = e11_codes_first)
icd_min_date_indexes <- apply(icd_dates, 1, function(x) {
    ifelse(
        all(is.na(x)),
        NA,
        which(x == min(x, na.rm = TRUE))
    )
})
icd_date_first <- sapply(1:length(icd_min_date_indexes), function(i) {
    if (is.na(icd_min_date_indexes[[i]])) {
        NA
    } else {
        icd_dates[i, icd_min_date_indexes[[i]]]
    }
}) %>% as.Date
icd_codes_first <- sapply(1:length(icd_min_date_indexes), function(i) {
    if (is.na(icd_min_date_indexes[[i]])) {
        NA
    } else {
        icd_codes[i, icd_min_date_indexes[[i]]]
    }
})

death_date_udis <- translations_present %>% filter(Description == "Date of death") %>% pull(UDI)

death_date_table <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = death_date_udis) %>% as.data.frame

death_date_table <- lapply(death_date_table, as.Date) %>% as.data.frame

baseline_date <- as.Date("2014-02-01")

death_before_baseline <- death_date_table[, 1] <= baseline_date

# Total n = 502,359 (including deaths)
# Deaths before baseline n = 8,946
# Deaths after baseline n = 35,553
incident <- icd_date_first > baseline_date # n = 23,623
prevalent <- icd_date_first <= baseline_date # n = 20,827


date_and_code_df <- data.frame(eid = eids, date = icd_date_first, code = icd_codes_first)
# saveRDS(date_and_code_df, "../data/t2d_date_and_code_df.rds")