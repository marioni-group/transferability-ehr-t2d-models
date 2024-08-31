library(data.table)
library(dplyr)
library(Matrix)

eid_index_range <- c(500001, 550000)
eid_index <- eid_index_range[[1]]:eid_index_range[[2]]


df <- fread("/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/omop_data/omop_condition_occurrence.txt.gz", 
            sep = "\t",
            select = c("eid", "condition_concept_id", "condition_start_date")) %>% as.data.frame

concept_table <- fread("/scratch/project_2007428/projects/prj_008_t2d_ehr/data/omop_concept.tsv",
                       sep = "\t")

df <- df %>% left_join(concept_table %>% select(concept_id, concept_code, concept_name, vocabulary_id), by = c(condition_concept_id = "concept_id"))

# Filter to SNOMED only
df <- df %>% filter(vocabulary_id == "SNOMED")
df <- df %>% select(-vocabulary_id)

print(paste0("Number of events before filtering to events before baseline: ", nrow(df)))

# Filter by occurrence before baseline date
df$condition_start_date <- as.Date(df$condition_start_date, format = "%d/%m/%Y")
df <- df %>% filter(condition_start_date < as.Date("2014-02-01"))

print(paste0("Number of events after filtering to events before baseline: ", nrow(df)))

# Extract unique concept_names and create dummy participant at the top with all codes to ensure all columns are present in final table
unique_concept_names <- df %>% pull(concept_name) %>% unique %>% sort
# Create dummy rows, one for each concept_name
dummy_rows <- lapply(unique_concept_names, function(x) {
    data.frame(eid = -1, 
               condition_concept_id = -1, 
               condition_start_date = "01/01/2013", 
               concept_code = -1, 
               concept_name = x)
})
dummy_row_df <- do.call(rbind, dummy_rows)

# Combine dummy rows with original
df <- rbind(dummy_row_df, df)

dummy_rows <- NULL
dummy_row_df <- NULL
gc()



term_colname <- "concept_name" # "condition_concept_id"

df[, term_colname] <- as.character(df[, term_colname])

# eid_table <- df %>% select(eid) %>% unique

eid_table <- fread("../../../data/base_data/ukbb_78537/phenotypes/ukbb_78537_base_phenotypes.csv", select = "eid")

# Check if upper index range is greater than number of eids in phenotype table and truncate if so.
eid_index_range[[2]] <- min(nrow(eid_table), eid_index_range[[2]])
eid_index <- eid_index_range[[1]]:eid_index_range[[2]]


# For testing: subset to 50,000 individuals
eid_table <- eid_table[eid_index, ,drop = FALSE]
# Add dummy
eid_table <- rbind(data.frame(eid = -1), eid_table)
df <- df %>% filter(eid %in% eid_table$eid) # Here the -1 ensures that the dummy participant is included

# First group by concept, then split into list of tbls
group_split_result <- df %>% 
    group_by(across(all_of(term_colname))) %>%
    select(eid, all_of(term_colname)) %>%
    group_split

terms <- sapply(group_split_result, function(term_df) {
    pull(term_df, term_colname)[[1]]
})

counts <- sapply(group_split_result, function(term_df) {
    nrow(term_df)
})
names(counts) <- terms

# For each concept, create a table of counts - one row per eid (participant) with the corresponding concept
# Then left_join each resulting count table with the original eid table for consistent row ids across the table
result <- lapply(group_split_result, function(term_df) {
    term <- pull(term_df, term_colname)[[1]]
    term_df <- term_df %>% count(eid) %>% rename('{{term}}' := n)
    eid_table %>% left_join(term_df, by = 'eid') %>% select(-eid)
})

# cbind all the resulting (single column) dataframes together
result <- bind_cols(result)
result <- bind_cols(list(eid_table, result)) %>% replace(is.na(.), 0)
result <- as(result, "sparseMatrix")
colnames(result) <- gsub('"', '', colnames(result))
saveRDS(result, paste0("../data/condition_frequency_table_eid_index_", eid_index_range[[1]], "-", eid_index_range[[2]], "_with_dummy.rds"))
