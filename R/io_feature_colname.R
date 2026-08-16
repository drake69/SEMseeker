#' Compose the column name of a computed feature (internal)
#'
#' AI-248. **The** single compositor of column names for the per-sample
#' statistics sibling, used by the producer that writes it and by the consumer
#' that references it. If the two sides built the string independently the
#' contract would break at the first divergence, so nothing else in the package
#' is allowed to paste these names together.
#'
#' The taxonomy has four axes, and the name carries them in order:
#'
#' \preformatted{
#' <SCOPE>_<MARKER>_<FIGURE>_<AGGREGATION>
#'
#' SCOPE        SAMPLE (whole sample) | <AREA> | <AREA>_<SUBAREA>
#' MARKER       SIGNAL | MUTATIONS | LESIONS | DELTA*
#' FIGURE       HYPER / HYPO for the instability markers,
#'              BETA / MVALUE (the scale) for SIGNAL
#' AGGREGATION  SUM | MEAN | MEDIAN | VARIANCE | IQR | MODELOW | MODEHIGH
#' }
#'
#' Until AI-248 the aggregation was implicit — one operator per marker, so
#' `MARKER`/`FIGURE` identified the value. With several aggregations over the
#' same scope it no longer does, hence the fourth axis.
#'
#' Axes left `NULL` are dropped, which is how scope-level properties are named:
#' `io_feature_colname("SAMPLE", aggregation = "N_PROBES")` gives
#' `SAMPLE_N_PROBES` — the number of positions is a property of the scope, not
#' an aggregation of a marker.
#'
#' @param scope scope name, as returned by [io_scope_name()]. Required.
#' @param marker,figure,aggregation the remaining axes; `NULL` to omit.
#' @return the column name, uppercase via [core_name_cleaning()].
#' @keywords internal
#' @noRd
io_feature_colname <- function(scope, marker = NULL, figure = NULL,
                               aggregation = NULL) {
  if (missing(scope) || is.null(scope) || !all(nzchar(scope)))
    stop("io_feature_colname(): 'scope' is required (use io_scope_name()).")

  parts <- list(scope, marker, figure, aggregation)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  parts <- lapply(parts, as.character)

  core_name_cleaning(do.call(paste, c(parts, list(sep = "_"))))
}
