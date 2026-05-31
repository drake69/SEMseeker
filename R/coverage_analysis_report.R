#' @title Probe coverage analysis report
#' @description Generate a coverage analysis report for Illumina methylation
#'   array data, summarising probe representation across genomic regions.
#'   WGBS data are not supported and will cause the function to stop with an
#'   informative error.
#' @param signal_data character. Path to the signal parquet file (e.g.
#'   \code{Data/SIGNAL_MEAN_PROBE_WHOLE.parquet}) or a data.frame already
#'   loaded into memory.
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend passed to
#'   \code{future} (default \code{"multicore"}).
#' @param ... Additional named arguments passed to \code{init_env()}.
#' @return Invisibly \code{NULL}. Coverage tables and charts are written to
#'   the result folder.
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' coverage_analysis_report(
#'   signal_data   = "~/semseeker_results/Data/SIGNAL_MEAN_PROBE_WHOLE.parquet",
#'   result_folder = "~/semseeker_results/"
#' )
#' }
#' @export
coverage_analysis_report <- function (signal_data, result_folder, maxResources = 90, parallel_strategy  = "multicore", ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, ...)
  # signal_data_path <- pivot_file_name_parquet("SIGNAL","MEAN","PROBE","WHOLE")
  # if (!file.exists(signal_data_path))
  # {
  #   log_event("ERROR:  ", format(Sys.time(), "%a %b %d %X %Y"), " Signal data is missing")
  #   close_env()
  #   stop()
  # }
  # signal_data <- polars::pl$read_parquet(signal_data_path)$to_data_frame()
  signal_data <- source_data_get(signal_data)
  ssEnv <- get_meth_tech(signal_data)
  if (ssEnv$tech=="WGBS")
  {
    log_event("ERROR:  ", format(Sys.time(), "%a %b %d %X %Y"), " WGBS data is not supported for coverage analysis")
    close_env()
    stop()
  }
  coverage_analysis(signal_data$AREA)
  close_env()
}
