#' Write the study sample sheet of the run (internal)
#'
#' AI-223: this function used to append the per-sample burden
#' (`MUTATIONS_HYPER`, `DELTAS_HYPO`, …) and `PROBES_COUNT` to the sample
#' sheet. Those columns were aggregated over the `AREA == "POSITION"` keys —
#' i.e. over every probe, with no genomic filter — which is exactly the
#' `SAMPLE` scope of the statistics sibling. They now live in
#' `SAMPLE_STATS_RESULT.csv` (see `sem_sample_stats_build()`), named
#' `SAMPLE_<MARKER>_<FIGURE>` / `SAMPLE_N_PROBES`, and the sample sheet stays
#' the description of the study. The two join on `Sample_ID`.
#'
#' @return Invisibly the path of the file written.
#' @keywords internal
#' @noRd
sem_study_summary_total <- function()
{
  ssEnv <- core_get_session_info()
  study_summary <- sem_study_summary_get(with_sample_stats = FALSE)

  if (is.null(study_summary) || nrow(study_summary) == 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No sample sheet available; SAMPLE_SHEET_RESULT not written.")
    return(invisible(NULL))
  }

  # io_file_path_build() uppercases via core_name_cleaning() — on-disk name is
  # SAMPLE_SHEET_RESULT.csv. Linux ext4 is case-sensitive; readers that
  # hard-code the path must use the uppercase form. See io_file_path_build().
  summary_file <- io_file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  utils::write.csv2(study_summary, summary_file, row.names = FALSE)

  invisible(summary_file)
}
