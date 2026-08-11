#' Reduce one sample's signal vector to its descriptors
#'
#' Generic descriptive statistics on a numeric vector: no SEM domain knowledge,
#' no I/O. Slice 1 of AI-223 calls it once per sample on the whole probe set
#' (scope `SAMPLE`); the region scopes reuse it unchanged on the subset of
#' probes of an area/subarea.
#'
#' The returned keys are built from [io_signal_stats()] rather than written out
#' here, so the computed statistics and the vocabulary declared to the consumer
#' cannot diverge. A statistic this function does not fill stays `NA`.
#'
#' Split at 0.5 for the two beta modes: methylation beta values are bimodal
#' with peaks near 0 (unmethylated) and near 1 (methylated), so the highest
#' density peak on each side of 0.5 is the natural estimate of each mode. On
#' M-values the distribution is unimodal and the split is meaningless, which is
#' why [io_signal_stats()] does not declare the mode columns there.
#'
#' @param values numeric vector, already filtered to finite values.
#' @param beta logical. TRUE when the signal is on the [0,1] beta scale.
#' @return named list of statistics, one element per [io_signal_stats()] entry.
#' @keywords internal
#' @noRd
util_signal_descriptors <- function(values, beta = TRUE) {

  wanted <- io_signal_stats(beta)
  out    <- stats::setNames(vector("list", length(wanted)), wanted)
  out[]  <- NA_real_

  out$N_PROBES <- length(values)

  if (length(values) == 0)
    return(out)

  out$MEDIAN   <- stats::median(values)
  out$MEAN     <- mean(values)
  out$VARIANCE <- stats::var(values)
  out$IQR      <- stats::IQR(values)

  # Off the beta scale the modes are not part of the vocabulary at all, so
  # there is nothing left to estimate.
  if (!all(c("MODE_LOW", "MODE_HIGH") %in% wanted))
    return(out)

  # density() needs at least two distinct points to estimate anything.
  if (length(unique(values)) < 2L)
    return(out)

  dens <- tryCatch(stats::density(values, from = 0, to = 1),
                   error = function(e) NULL)
  if (is.null(dens))
    return(out)

  low <- dens$x < 0.5
  if (any(low))
    out$MODE_LOW <- dens$x[low][which.max(dens$y[low])]
  if (any(!low))
    out$MODE_HIGH <- dens$x[!low][which.max(dens$y[!low])]

  out
}
