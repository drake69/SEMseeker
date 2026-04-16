#' Core SEMseeker pipeline (internal)
#'
#' Internal entry point used by the public \code{\link{semseeker}} dispatcher.
#' Accepts already-normalised signal data (matrix or data frame with probe-ID
#' rownames) and runs the full SEM analysis. Users should call
#' \code{\link{semseeker}} instead, which handles input normalisation,
#' M-value conversion, and tech/genome_build validation.
#'
#' @param sample_sheet Data frame (or list of data frames) with a \code{Sample_ID} column.
#' @param signal_data Methylation matrix or data frame with probe-ID rownames
#'   (or list of such objects, one per batch).
#' @param result_folder Output directory.
#' @param ... Additional arguments forwarded to \code{init_env()}.
#'
#' @return Invisibly NULL; writes output files to \code{result_folder}.
#' @keywords internal
#' @importFrom doRNG %dorng%
semseeker_core <- function(sample_sheet,
  signal_data,
  result_folder,
  ... ) {

  init_env( result_folder= result_folder, ...)

  ssEnv <- get_session_info()
  log_event("BANNER:", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will search MARKERS for project \n in ", ssEnv$result_folderData)

  # check if the input is a list of data frames
  if(!is.list(sample_sheet) | is.data.frame(sample_sheet))
    sample_sheet <- list(sample_sheet)
  if(!is.list(signal_data) | is.data.frame(signal_data))
    signal_data <- list(signal_data)


  batch_id <- 1
  ssEnv$batch_count <- length(sample_sheet)
  ssEnv <- update_session_info(ssEnv)

  # C-06: write session provenance metadata (genome_build, tech, version, …)
  total_sample_n <- sum(sapply(sample_sheet, nrow))
  session_metadata_write(result_folder, sample_n = total_sample_n)

  for(batch_id in seq_along(sample_sheet))
  {
    start_time <- Sys.time()
    ssEnv$running_batch_id <- batch_id
    ssEnv <- update_session_info(ssEnv)
    sample_sheet_local <- source_data_get(sample_sheet[[batch_id]])
    sample_sheet_local$Sample_ID <- name_cleaning(sample_sheet_local$Sample_ID)
    utils::write.csv2(sample_sheet_local, file = file_path_build(ssEnv$result_folderData, paste0(batch_id,"_sample_sheet_original"),"csv",FALSE))
    analyze_batch(source_data_get(signal_data[[batch_id]]), sample_sheet_local)
    create_position_pivots(sample_sheet_local,ssEnv$keys_markers_figures)
    log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), "Batch Executed in:", difftime(time1 = Sys.time(), time2= start_time,units = "mins") , " minutes.")
  }

  deltaX_get()
  study_summary_total()
  annotate_position_pivots()
  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " Saving Sample Sheet with Results! ")

  close_env()
}
