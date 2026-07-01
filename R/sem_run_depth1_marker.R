#' Run depth=1 (sample-level) association for one marker
#'
#' Extracted from association_analysis() (was inline at lines 218-290).
#' Reads per-sample counts already present in study_summary columns,
#' optionally resumes from a partial results CSV, then applies the
#' stat model per key (one row per genomic region key).
#'
#' @param prep list returned by sem_prepare_study_for_analysis().
#' @param keys data.frame of keys for this marker (subset of
#'   ssEnv$keys_markers_figures).
#' @param family_test character.
#' @param fileNameResults character. Path of the output CSV.
#' @param filter_p_value logical.
#' @param ssEnv list. Session environment from core_get_session_info().
#' @param ... forwarded to apply_stat_model().
#' @return list(results = data.frame, processed_items = integer).
#'   Side effect: writes the CSV via association_analysis_save_results().
#' @keywords internal
sem_run_depth1_marker <- function(prep, keys, family_test, fileNameResults,
                               filter_p_value, ssEnv, ...) {
  results <- data.frame()
  processed_items <- 0L

  cols <- keys$COMBINED
  if (sum(cols %in% colnames(prep$study_summary)) == 0)
    return(list(results = results, processed_items = processed_items))

  cols          <- cols[cols %in% colnames(prep$study_summary)]
  keys$AREA     <- "SAMPLE_GROUP"
  keys$SUBAREA  <- "SAMPLE"

  has_covariates <- !is.null(prep$covariates) && length(prep$covariates) != 0
  if (has_covariates) {
    study_summary_local <- prep$study_summary[,
      c(prep$independent_variable, prep$covariates, cols, "Sample_Group")]
  } else {
    study_summary_local <- prep$study_summary
  }

  # resume: drop keys already present in the CSV (DEPTH==1 only)
  file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size > 3
  old_results <- data.frame()
  if (file_good) {
    old_results <- unique(utils::read.csv2(fileNameResults, header = TRUE))
    old_filtered <- old_results[old_results$DEPTH == 1, ]
    done_ids <- unlist(apply(unique(old_filtered[, c("MARKER", "FIGURE", "AREA", "SUBAREA")]), 1,
      function(x) paste(x, collapse = "_", sep = "")))
    todo_ids <- unlist(apply(keys[, c("MARKER", "FIGURE", "AREA", "SUBAREA")], 1,
      function(x) paste(x, collapse = "_", sep = "")))
    keys <- keys[!(todo_ids %in% done_ids), ]
  }

  if (nrow(keys) > 0) {
    for (j in seq_len(nrow(keys))) {
      key <- keys[j, ]
      key$FIGURE <- as.character(key$FIGURE)
      key$MARKER <- as.character(key$MARKER)
      g_start <- 2 + length(prep$covariates)
      column_selectors <- c(prep$independent_variable, prep$covariates, key$COMBINED)
      column_selectors <- column_selectors[column_selectors != ""]
      processed_items <- processed_items + ncol(study_summary_local) - g_start
      if (any(is.na(study_summary_local[, column_selectors]))) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
          " Missing values in the data frame!")
        study_summary_local <- study_summary_local[
          complete.cases(study_summary_local[, column_selectors]), ]
      }
      result_temp <- apply_stat_model(
        tempDataFrame  = study_summary_local[, column_selectors],
        g_start        = g_start,
        family_test    = family_test,
        covariates     = prep$covariates,
        key            = key,
        transformation_y = prep$transformation_y,
        dototal        = FALSE,
        session_folder = ssEnv$session_folder,
        prep$independent_variable,
        prep$depth_analysis,
        prep$inference_detail$samples_sql_condition,
        inference_detail = prep$inference_detail,
        ...)
      results <- plyr::rbind.fill(results, result_temp)
    }
  }

  if (!is.null(dim(results)) && nrow(results) > 0 && "PVALUE_ADJ" %in% colnames(results))
    results <- results[order(results$PVALUE_ADJ), ]

  if (nrow(old_results) > 0) {
    results <- plyr::rbind.fill(results, old_results)
  }

  results[results == ""] <- NA
  if (ncol(results) > 0)
    results <- results[, colSums(is.na(results)) < nrow(results), drop = FALSE]
  results <- results[, !grepl("SAMPLES_SQL_CONDITION", colnames(results)), drop = FALSE]
  association_analysis_save_results(results, fileNameResults, family_test, filter_p_value)

  list(results = results, processed_items = processed_items)
}
