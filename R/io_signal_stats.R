#' Statistics emitted for a signal scope
#'
#' Single source of truth for which statistics a scope carries: the producer
#' composes the column names from this vector via [io_stat_colname()], and
#' [util_signal_descriptors()] builds its output on it, so the computed keys and
#' the declared vocabulary cannot drift apart.
#'
#' MODE_LOW / MODE_HIGH only make sense on a bounded, bimodal scale (beta in
#' [0,1], the two peaks around 0 and 1). On M-values the distribution is
#' roughly gaussian and the two-mode split carries no meaning, so the columns
#' are omitted rather than filled with a meaningless number.
#'
#' @param beta logical. TRUE when the signal is on the [0,1] beta scale.
#' @return character vector of statistic names.
#' @keywords internal
#' @noRd
io_signal_stats <- function(beta = TRUE) {
  base <- c("MEDIAN", "MEAN", "VARIANCE", "IQR", "N_PROBES")
  if (isTRUE(beta)) c(base, "MODE_LOW", "MODE_HIGH") else base
}
