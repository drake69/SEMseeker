#' Calculate Stochastic Epigenetic Mutations from a Methylation Dataset
#'
#' The `semseeker` function processes a methylation dataset to identify stochastic epigenetic mutations and generates output reports. This involves working with sample sheets and signal data to produce pivot tables and bedgraph files.
#'
#' @param sample_sheet A data frame containing at least a column named `Sample_ID` to identify samples. This can be a single data frame or a list of data frames.
#' @param signal_data A matrix of methylation data. This can be a single data frame or a list of data frames.
#' @param result_folder A string specifying the directory where the results will be saved.
#' @param ... Additional arguments for further processing options, including:
#'   - `parallel_strategy`: Strategy for parallel execution. Possible values are `none`, `multisession`, `sequential`, `multicore`, and `cluster`.
#'   - `maxResources`: Percentage of available cores to be used, default is 90 percent, rounded to the lowest integer.
#'   - `signal_intrasample`: A logical value indicating whether the signal data is intrasample. Default is `FALSE`.
#'   - `sliding_window_size`: An integer specifying the size of the sliding window. Default is 11.
#'   - `alpha`: A numeric value specifying the alpha threshold for the analysis. Default is 0.05.
#'   - `showprogress`: A logical value indicating whether to show a progress bar. Default is `TRUE`.
#'   - `iqrTimes`: A numeric value specifying the interquartile range multiplier to identify aberration in the data. Default is 3.
#'   - `sex_chromosome_remove`: A logical value indicating whether to remove the sex chromosomes. Default is `TRUE`.
#'   - `plot_format`: A string specifying the plot format. Default is "png".
#'   - `verbosity`: A numeric value specifying the verbosity level. Default is 0.
#'   - `marker`: A vector of string specifying the marker to be used for the analysis. Default are "MUTATIONS", "LESIONS","DELTARQ","DELTAQ","DELTAS","DELTAR","MEAN","DELTARP","DELTAP"
#'   - `areas`: A vector of string specifying the areas to be used for the analysis. Default are "GENE","CHR","ISLAND","PROBE"
#'
#' @return The function writes multiple files to the specified `result_folder`, including the processed sample sheet and result files in CSV format, along with pivot tables and bedgraph files.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' sample_sheet <- data.frame(
#'   Sample_ID    = c("S1", "S2", "S3"),
#'   Sample_Name  = c("Case1", "Case2", "Ctrl1"),
#'   Sample_Group = c("Case", "Case", "Control")
#' )
#' semseeker(
#'   sample_sheet  = sample_sheet,
#'   signal_data   = beta_matrix,
#'   result_folder = tempdir()
#' )
#' }
#' @export
#' @importFrom doRNG %dorng%
#'
semseeker <- function(sample_sheet,
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
    signal_intrasample <- TRUE
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
