#' Final per-job save with TRANSFORMATION_X annotation
#'
#' Extracted from association_analysis() (was inline at lines 436-444).
#' Adds the TRANSFORMATION_X column to the last marker's results and
#' performs the final CSV save; then emits the JOURNAL closing entry.
#'
#' @param results data.frame from the last marker iteration (may be empty).
#' @param inference_detail single-row data.frame.
#' @param family_test character.
#' @param filter_p_value logical.
#' @param fileNameResults character. Path of the output CSV.
#' @param start_time POSIXct. Job start.
#' @param processed_items integer. Total items processed in this job.
#' @return Invisibly NULL.
#' @keywords internal
finalize_job_results <- function(results, inference_detail, family_test,
                                  filter_p_value, fileNameResults,
                                  start_time, processed_items) {
  if (!is.null(results) && nrow(results) != 0) {
    results$TRANSFORMATION_X <- inference_detail$transformation_x
    association_analysis_save_results(results, fileNameResults, family_test, filter_p_value)
  }
  total_time <- difftime(Sys.time(), start_time, units = "mins")
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
    " Finished processing association analysis for ", processed_items,
    " items in ", total_time, " minutes.")
  log_event("JOURNAL:", format(Sys.time(), "%a %b %d %X %Y"),
    " Association Analysis finished in ", total_time,
    " minutes. \n ####################################################")
  invisible(NULL)
}
