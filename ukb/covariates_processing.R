library(dplyr)
library(data.table)
library(zoo)

phenotypes_data_path <- "../../../data/base_data/ukbb_78537/phenotypes/"

translations <- read.csv(paste0(phenotypes_data_path, "ukbb_78537_translations.tsv"), sep = "\t")

base_phenotypes_head <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), nrows = 5)

base_phenotypes_colnames <- colnames(base_phenotypes_head)

# Filter translations to those corresponding to columns that are present in base_phenotypes
translations_present <- translations %>% filter(UDI %in% base_phenotypes_colnames)

udis <- c("34-0.0", "52-0.0", "31-0.0")
covariate_names <- c("yob", "mob", "sex")

eids <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = "eid")

covariates <- fread(paste0(phenotypes_data_path, "ukbb_78537_base_phenotypes.csv"), select = udis) %>% as.data.frame

colnames(covariates) <- covariate_names

covariates$yobmob <- paste0(as.character(covariates$yob), "-", as.character(covariates$mob))
covariates$yobmob <- as.yearmon(covariates$yobmob)

covariates$age <- as.yearmon("2014-02") - covariates$yobmob

# Sex is female=1, male=0 in FinnGen. Flip the UKB coding to match
covariates$sex <- 1 - covariates$sex

covariates$eid <- eids$eid
covariates <- covariates %>% select(eid, age, sex)

saveRDS(covariates, "../data/covariates.rds")
