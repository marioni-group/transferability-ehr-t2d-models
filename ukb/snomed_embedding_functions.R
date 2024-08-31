library(dplyr)

load_embedding_table <- function() {
  embedding_table <- unique(read.csv("../data/embeddings.csv"))
  row.names(embedding_table) <- embedding_table[, "X"]
  embedding_table[, "X"] <- NULL
  embedding_table
}

load_icd_embedding_table <- function() {
  icd_embeddings <- read.csv("../data/all_icd10_embeddings_text-embedding-ada-002.csv")

  embedding_names <- icd_embeddings[, 1]

  icd_embeddings <- icd_embeddings[!duplicated(embedding_names), ]

  row.names(icd_embeddings) <- icd_embeddings[, 1]
  icd_embeddings <- icd_embeddings[, -1]
  colnames(icd_embeddings) <- sapply(0:(ncol(icd_embeddings) - 1), function(x) {paste0("embdim_", x)})
  icd_embeddings
}

remove_disorder_or_finding_suffix <- function(terms) {
  terms %>% sapply(function(x) {gsub(" \\(disorder\\)", "", x)}) %>% sapply(function(x) {gsub(" \\(finding\\)", "", x)})
}

# Assumes row.names of embedding_table are SNOMED terms
get_embedding <- function(term, embedding_table) {
  embedding_table[term, ]
}


find_nearest_present_term <- function(nn2_result, embedding_table, check_presence_terms) {
  get_first_present_term_concept <- function(embedding_table_indexes) {
    terms <- embedding_table[embedding_table_indexes, ] %>% row.names
    concepts <- terms %>% sapply(term_to_concept)
    presence <- concepts %in% check_presence_terms
    first_presence_index <- which(presence)[1]
    result_term <- ifelse(is.na(first_presence_index), NA, terms[[first_presence_index]])
    result_concept <- ifelse(is.na(first_presence_index), NA, concepts[[first_presence_index]])
    list(term = result_term, concept = result_concept)
  }

  apply(nn2_result$nn.idx, 1, get_first_present_term_concept)
}

