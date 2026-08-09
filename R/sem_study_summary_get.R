#' Read the study sample sheet, optionally joined with the statistics sibling
#'
#' @param sql_sample_selection filter applied via [assoc_filter_sql()].
#' @param with_sample_stats logical. When `TRUE` (default) the per-sample
#'   statistics produced by `sem_sample_stats_build()` are LEFT JOINed on
#'   `Sample_ID`, so consumers see `SAMPLE_<MARKER>_<FIGURE>` and
#'   `SAMPLE_<STAT>` as ordinary columns (AI-223). The producers of the two
#'   files call it with `FALSE` to avoid feeding their own output back in.
#' @keywords internal
#' @noRd
sem_study_summary_get <- function(sql_sample_selection="", with_sample_stats = TRUE)
{
  ssEnv <- core_get_session_info()
  # io_file_path_build() uppercases via core_name_cleaning() — the on-disk file is
  # SAMPLE_SHEET_RESULT.csv, written by sem_study_summary_total(). Going through
  # io_file_path_build() here guarantees the same path is resolved on every OS
  # regardless of file-system case sensitivity.
  summary_file <- io_file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  if(!file.exists(summary_file))
    summary_file <- io_file_path_build( ssEnv$result_folderData, "1_SAMPLE_SHEET_ORIGINAL","csv")

  if(!file.exists(summary_file))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No sample sheet found in ", ssEnv$result_folderData)
    return()
  }

  study_summary <-   utils::read.csv2(summary_file)
  study_summary <- assoc_filter_sql(sql_sample_selection, study_summary)

  reference <- study_summary[study_summary$Sample_Group == "Reference",]
  no_reference <- study_summary[study_summary$Sample_Group != "Reference",]
  reference <- reference[!(reference$Sample_ID %in% no_reference$Sample_ID), ]

  study_summary <- rbind(no_reference, reference)

  if(nrow(reference)==0)
    study_summary <- no_reference

  # AI-223: the burden and the signal descriptors live in the sibling; join
  # them here so every consumer keeps seeing one table.
  if (isTRUE(with_sample_stats))
    study_summary <- .sem_study_summary_join_stats(study_summary)

  return(study_summary)
}

#' LEFT JOIN of the per-sample statistics sibling onto the sample sheet
#'
#' Silent no-op when the sibling has not been produced (e.g. a session that
#' only ran the exploratory step, or a result folder from an older version):
#' the sample sheet is still perfectly usable on its own, callers that need the
#' burden check for their columns.
#'
#' @keywords internal
#' @noRd
.sem_study_summary_join_stats <- function(study_summary) {

  if (is.null(study_summary) || !("Sample_ID" %in% colnames(study_summary)))
    return(study_summary)

  ssEnv <- core_get_session_info()
  stats_file <- io_file_path_build(ssEnv$result_folderData, "sample_stats_result", "csv")
  if (!file.exists(stats_file))
    return(study_summary)

  stats <- tryCatch(utils::read.csv2(stats_file, stringsAsFactors = FALSE),
                    error = function(e) NULL)
  if (is.null(stats) || !("Sample_ID" %in% colnames(stats)) || nrow(stats) == 0)
    return(study_summary)

  # Never let the sibling shadow a column of the sample sheet.
  duplicated_cols <- setdiff(intersect(colnames(stats), colnames(study_summary)), "Sample_ID")
  if (length(duplicated_cols) > 0)
    stats <- stats[, !(colnames(stats) %in% duplicated_cols), drop = FALSE]

  merge(study_summary, stats, by = "Sample_ID", all.x = TRUE)
}
