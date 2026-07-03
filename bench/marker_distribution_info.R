#' Marker distribution charts and theoretical fit
#'
#' Generates per-marker distribution diagnostics for a completed SEMseeker
#' run: per-marker quantisation metrics + histogram/density plots with
#' optional theoretical-distribution fits (Poisson, Negative Binomial,
#' Gaussian, ...). Output charts are written under
#' \code{<result_folder>/Chart/Distributions/}.
#'
#' @param result_folder character. Path to the SEMseeker result folder of
#'   a completed run (must contain \code{Data/Pivots/}).
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend passed to
#'   \code{future} (default \code{"multisession"}).
#' @param ... Additional named arguments passed to \code{init_env()}.
#'
#' @return Invisibly \code{NULL}. PNG charts and a quantisation-metrics
#'   table are written to disk as a side effect.
#'
#' @keywords internal
#' @noRd
marker_distribution_info <- function(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = FALSE, ...)

  # 
  # 

  marker_quantization_metric(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
  marker_fit_to_theoretical_distribution()

  close_env()

}
