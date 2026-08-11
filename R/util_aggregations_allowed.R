#' Aggregations admissible for a marker/figure pair (internal)
#'
#' AI-248. Single source of truth of the AGGREGATION axis: which ways of
#' reducing a set of positions to one number are admissible for a given
#' `(MARKER, FIGURE)`. The producer iterates over this to know what to compute,
#' the consumer validates against it, and [io_feature_colname()] names the
#' result — so the computed keys and the declared vocabulary cannot drift apart.
#'
#' The axis is **permissive by design**: every aggregation applies to every
#' marker. Whether the median of hypermethylated mutations *means* something is
#' a question for the biologist reading the table; making it expressible and
#' nameable in a uniform way is the package's job. This generalises the old
#' binary `DISCRETE` flag, which could only say "sum" or "mean".
#'
#' One restriction, and it is not arbitrary: `MODE_LOW` / `MODE_HIGH` estimate
#' the two peaks of a bimodal distribution on a bounded scale, which only exists
#' for `SIGNAL` on the `BETA` scale. On M-values the distribution is unimodal and
#' the split carries no meaning.
#'
#' @param marker,figure the pair. `figure` matters only for `SIGNAL`, where it
#'   carries the scale (see [io_signal_figure()]).
#' @param discrete the `DISCRETE` flag of the key. It no longer decides the
#'   operator — that is the point of this axis — but it still says which
#'   aggregation is produced **by default** for a marker that is not `SIGNAL`:
#'   a count is summed, a continuous deviation is averaged. Preserving that
#'   distinction is what keeps the migration from silently changing the numbers
#'   of the DELTA family.
#' @param default when `TRUE` (default) return only what is produced without
#'   being asked. When `FALSE` return everything the taxonomy admits.
#' @return character vector of aggregation names.
#' @keywords internal
#' @noRd
util_aggregations_allowed <- function(marker, figure = NULL, discrete = TRUE,
                                      default = TRUE) {

  marker <- toupper(as.character(marker))
  figure <- if (is.null(figure)) "" else toupper(as.character(figure))

  is_signal <- identical(marker, "SIGNAL")
  bounded   <- is_signal && identical(figure, "BETA")

  every <- c("SUM", "MEAN", "MEDIAN", "VARIANCE", "IQR")
  if (bounded)
    every <- c(every, "MODE_LOW", "MODE_HIGH")

  if (!isTRUE(default))
    return(every)

  # Produced by default: the operator that carries the meaning of the marker.
  # For the raw signal there is no privileged operator, so the whole descriptor
  # set is produced. For everything else the default is the operator the package
  # has always applied — a count of epimutations is a burden and the burden is
  # its sum, a continuous deviation is averaged.
  if (is_signal)
    return(setdiff(every, "SUM"))

  if (isTRUE(discrete)) "SUM" else "MEAN"
}
