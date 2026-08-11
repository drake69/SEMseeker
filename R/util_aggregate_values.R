#' Reduce a set of values to one number, by aggregation name (internal)
#'
#' AI-248. The other half of the AGGREGATION axis: [util_aggregations_allowed()]
#' declares which names are legal, this one is the only place that says what each
#' name *computes*. Keeping the two together in one vocabulary is what stops a
#' column called `MEDIAN` from holding a mean.
#'
#' `MODE_LOW` / `MODE_HIGH` are the only distribution-shape operators: they
#' estimate the highest density peak on each side of 0.5, so they need the whole
#' distribution rather than a streaming reduce, and they are admissible only on
#' the bounded scale. They delegate to [util_signal_descriptors()], which already
#' implements the estimate.
#'
#' @param values numeric vector. Non-finite values are dropped.
#' @param aggregation one name from [util_aggregations_allowed()], or
#'   `"N_PROBES"` for the count of usable positions.
#' @return a single numeric, `NA_real_` when the values do not support it.
#' @keywords internal
#' @noRd
util_aggregate_values <- function(values, aggregation) {

  aggregation <- toupper(as.character(aggregation))
  values <- as.numeric(values)
  values <- values[is.finite(values)]

  if (identical(aggregation, "N_PROBES"))
    return(length(values))

  if (length(values) == 0L)
    return(NA_real_)

  switch(
    aggregation,
    SUM       = sum(values),
    MEAN      = mean(values),
    MEDIAN    = stats::median(values),
    VARIANCE  = stats::var(values),
    IQR       = stats::IQR(values),
    MODE_LOW  = util_signal_descriptors(values, beta = TRUE)$MODE_LOW,
    MODE_HIGH = util_signal_descriptors(values, beta = TRUE)$MODE_HIGH,
    stop("util_aggregate_values(): unknown aggregation '", aggregation,
         "'. Legal names come from util_aggregations_allowed().")
  )
}
