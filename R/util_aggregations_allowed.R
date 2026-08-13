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

  # Two classes, and the class is what decides. SIGNAL, DELTAS and DELTAR are
  # all continuous — a value per position — and take the same operators; the
  # instability markers are counts, 0 or 1 per position.
  # SIGNAL is continuous by definition, whatever the caller passes: the raw
  # signal is never a count.
  if (identical(marker, "SIGNAL"))
    discrete <- FALSE
  bounded <- identical(marker, "SIGNAL") && identical(figure, "BETA")

  if (isTRUE(discrete)) {
    # SUM is the burden. MEAN of a 0/1 vector is the same burden divided by the
    # positions of the scope, i.e. the DENSITY of epimutations — the only form
    # comparable across regions of different size, which a raw burden is not.
    # MEDIAN, VARIANCE and IQR of a 0/1 vector are degenerate, so they are not
    # admissible: naming them would only add noise to the file.
    admissible <- c("SUM", "MEAN")
    return(if (isTRUE(default)) "SUM" else admissible)
  }

  admissible <- c("SUM", "MEAN", "MEDIAN", "VARIANCE", "IQR")
  # The two modes look for the highest density peak on EACH SIDE OF 0.5, which
  # presupposes a scale bounded in [0,1] and bimodal. DELTAS and DELTAR are
  # deviations centred on zero: that split says nothing about them.
  if (bounded)
    admissible <- c(admissible, "MODE_LOW", "MODE_HIGH")

  if (!isTRUE(default))
    return(admissible)

  # Every continuous marker is produced with the whole descriptor set: no
  # operator is privileged for a continuous value. This ADDS columns without
  # changing any existing number — the historical mean keeps its value and
  # merely gains the _MEAN suffix.
  setdiff(admissible, "SUM")
}
