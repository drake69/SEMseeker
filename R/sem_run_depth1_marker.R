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
#' @param ... forwarded to assoc_apply_stat_model().
#' @return list(results = data.frame, processed_items = integer).
#'   Side effect: writes the CSV via assoc_analysis_save_results().
#' @keywords internal
sem_run_depth1_marker <- function(prep, keys, family_test, fileNameResults,
                               filter_p_value, ssEnv, ...) {
  results <- data.frame()
  processed_items <- 0L

  # AI-223: the sample-level burden lives in the statistics sibling, joined
  # onto study_summary by sem_study_summary_get(). Its column names are
  # composed with the shared helper — never hard-coded — and renamed back to
  # <MARKER>_<FIGURE> below so the inference output (AREA_OF_TEST) keeps its
  # historical values.
  # AI-223 slice 2a: the same test can run on a region scope of the sibling
  # (e.g. GENE_TSS1500 = burden restricted to the probes of that region class).
  # It is still one value per sample, hence still DEPTH=1; the scope is what
  # the AREA column carries.
  scopes <- .sem_depth1_scopes_get(prep, ssEnv)

  # resume: drop keys already present in the CSV (DEPTH==1 only). Read once,
  # the scope is part of the identity through AREA.
  file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size > 3
  old_results <- data.frame()
  done_ids <- character(0)
  if (file_good) {
    old_results <- unique(utils::read.csv2(fileNameResults, header = TRUE))
    old_filtered <- old_results[old_results$DEPTH == 1, ]
    done_ids <- unlist(apply(unique(old_filtered[, c("MARKER", "FIGURE", "AREA", "SUBAREA")]), 1,
      function(x) paste(x, collapse = "_", sep = "")))
  }

  any_scope_ran <- FALSE
  for (scope in scopes) {
    scope_keys  <- keys
    cols        <- scope_keys$COMBINED
    burden_cols <- io_burden_colname(scope, scope_keys$MARKER, scope_keys$FIGURE)
    present     <- burden_cols %in% colnames(prep$study_summary)
    if (sum(present) == 0) {
      if (!identical(scope, "SAMPLE"))
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
          " Scope ", scope, " carries no burden column for this marker; skipped.")
      next
    }

    any_scope_ran <- TRUE
    scope_keys  <- scope_keys[present, , drop = FALSE]
    cols        <- cols[present]
    burden_cols <- burden_cols[present]
    # The renaming is per scope and starts from the untouched study_summary:
    # two scopes would otherwise fight over the same <MARKER>_<FIGURE> name.
    study_summary_scope <- .sem_burden_cols_rename(prep$study_summary, burden_cols, cols)
    # SUBAREA == "SAMPLE" is what makes assoc_analysis_save_results() stamp
    # DEPTH = 1; AREA says which scope was aggregated.
    scope_keys$AREA    <- if (identical(scope, "SAMPLE")) "SAMPLE_GROUP" else scope
    scope_keys$SUBAREA <- "SAMPLE"

    has_covariates <- !is.null(prep$covariates) && length(prep$covariates) != 0
    if (has_covariates) {
      study_summary_local <- study_summary_scope[,
        c(prep$independent_variable, prep$covariates, cols, "Sample_Group")]
    } else {
      study_summary_local <- study_summary_scope
    }

    if (length(done_ids) > 0) {
      todo_ids <- unlist(apply(scope_keys[, c("MARKER", "FIGURE", "AREA", "SUBAREA")], 1,
        function(x) paste(x, collapse = "_", sep = "")))
      scope_keys <- scope_keys[!(todo_ids %in% done_ids), ]
    }

    if (nrow(scope_keys) == 0)
      next

    for (j in seq_len(nrow(scope_keys))) {
      key <- scope_keys[j, ]
      key$FIGURE <- as.character(key$FIGURE)
      key$MARKER <- as.character(key$MARKER)
      g_start <- 2 + length(prep$covariates)
      column_selectors <- c(prep$independent_variable, prep$covariates, key$COMBINED)
      column_selectors <- column_selectors[column_selectors != ""]
      # One key at depth=1 means exactly one test: a single burden column is
      # fitted against the independent variable. The counter used to take
      # ncol(study_summary_local), i.e. the width of the WHOLE sample sheet,
      # which counted phenotype columns as if they had been tested (and, since
      # AI-223 joins the statistics sibling, the signal descriptors too).
      # Depth>1 counts the tested columns of its batch frame; this is the
      # depth=1 equivalent.
      processed_items <- processed_items + 1L
      if (any(is.na(study_summary_local[, column_selectors]))) {
        core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
          " Missing values in the data frame!")
        study_summary_local <- study_summary_local[
          complete.cases(study_summary_local[, column_selectors]), ]
      }
      result_temp <- assoc_apply_stat_model(
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

  # Nothing of this marker is present in the sibling (e.g. a result folder
  # produced before AI-223): leave the CSV exactly as it was found.
  if (!any_scope_ran)
    return(list(results = data.frame(), processed_items = processed_items))

  if (!is.null(dim(results)) && nrow(results) > 0 && "PVALUE_ADJ" %in% colnames(results))
    results <- results[order(results$PVALUE_ADJ), ]

  if (nrow(old_results) > 0) {
    results <- plyr::rbind.fill(results, old_results)
  }

  results[results == ""] <- NA
  if (ncol(results) > 0)
    results <- results[, colSums(is.na(results)) < nrow(results), drop = FALSE]
  results <- results[, !grepl("SAMPLES_SQL_CONDITION", colnames(results)), drop = FALSE]
  assoc_analysis_save_results(results, fileNameResults, family_test, filter_p_value)

  list(results = results, processed_items = processed_items)
}

#' Scopes to test at depth=1, validated against the statistics sibling
#'
#' AI-223 slice 2a. `inference_details$scopes` names them the way the sibling
#' names its columns (`SAMPLE`, `GENE_TSS1500`, …), several separated by `"+"`
#' like covariates. Default `SAMPLE`, which is the historical behaviour.
#'
#' A requested scope whose columns are absent from the sibling is an error, not
#' a skip: the run would otherwise produce a complete-looking inference CSV
#' that simply never tested what was asked. The producer side of the contract
#' is `semseeker(sample_stats_scopes = ...)`.
#'
#' @keywords internal
#' @noRd
.sem_depth1_scopes_get <- function(prep, ssEnv) {

  requested <- prep$inference_detail$scopes
  if (is.null(requested) || length(requested) == 0 || all(is.na(requested)))
    return("SAMPLE")
  requested <- util_split_and_clean(requested)
  requested <- core_name_cleaning(requested[nzchar(requested)])
  if (length(requested) == 0)
    return("SAMPLE")

  # Every marker/figure the run knows about — the sibling carries a column per
  # pair and per scope, so one hit is enough to say the scope was produced.
  all_keys <- ssEnv$keys_markers_figures
  known_cols <- colnames(prep$study_summary)
  for (scope in requested) {
    candidates <- io_burden_colname(scope, all_keys$MARKER, all_keys$FIGURE)
    if (!any(candidates %in% known_cols))
      stop("inference_details$scopes: scope '", scope, "' has no column in ",
           "SAMPLE_STATS_RESULT. Produce it with ",
           "semseeker(sample_stats_scopes = c(\"SAMPLE\", \"", scope, "\")) ",
           "and rerun the analysis.")
  }
  unique(requested)
}

#' Rename the SAMPLE-scope burden columns back to <MARKER>_<FIGURE>
#'
#' AI-223. The statistics sibling stores the sample-level burden as
#' `SAMPLE_<MARKER>_<FIGURE>`; the models and the inference output have always
#' referred to it as `<MARKER>_<FIGURE>` (it becomes AREA_OF_TEST in the
#' results). Renaming at the boundary keeps the artefact naming explicit about
#' its scope without changing what the results look like.
#'
#' @keywords internal
#' @noRd
.sem_burden_cols_rename <- function(study_summary, from, to) {
  if (is.null(study_summary) || length(from) == 0)
    return(study_summary)
  # drop any stale column that would collide with the target name
  collisions <- intersect(to, colnames(study_summary))
  if (length(collisions) > 0)
    study_summary <- study_summary[, !(colnames(study_summary) %in% collisions), drop = FALSE]
  idx <- match(from, colnames(study_summary))
  keep <- !is.na(idx)
  colnames(study_summary)[idx[keep]] <- to[keep]
  study_summary
}
