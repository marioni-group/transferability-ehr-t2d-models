library(SurvMetrics)
library(pec)

# This function assumes a result object with a summaries attribute - a list of summary objects, one for each model, indexed by summary_names
get_c_index <- function(result, summary_names = c("summary_test_null_model",
                                                  "summary_test_ehr_model",
                                                  "summary_test_prs_model",
                                                  "summary_test_full_model",
                                                  "summary_test_interaction_model"),
                        model_names = c("Null",
                                        "EHR",
                                        "PRS",
                                        "Full",
                                        "Interaction"),
                        new_colnames = c("C-Index", "SE")) {
  c_indexes <- lapply(summary_names, function(summary_name) {
    result[["summaries"]][[summary_name]]$concordance
  })
  c_indexes <- do.call(rbind, c_indexes)
  row.names(c_indexes) <- model_names
  colnames(c_indexes) <- new_colnames
  c_indexes
}

calculate_brier_scores <- function(result, result_model_names = c("null_model",
                                                                  "ehr_model",
                                                                  "prs_model",
                                                                  "full_model",
                                                                  "interaction_model"),
                                   model_names = c("Null",
                                                   "EHR",
                                                   "PRS",
                                                   "Full",
                                                   "Interaction"),
                                   ts = c(1:5)) {
  models <- lapply(result_model_names, function(result_model_name) {
    result[["cox_models"]][[result_model_name]]
  })
  brier_scores <- lapply(models, function(model) {
    sapply(ts, function(time_point) {
      survival_probability <- predictSurvProb(model, )
      Brier()
    })
  })
}