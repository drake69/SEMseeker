#' Every aggregation name the taxonomy knows (internal)
#'
#' AI-248. The full vocabulary of the AGGREGATION axis, independent of any
#' marker: the union of what [util_aggregations_allowed()] admits over both
#' classes of marker, counts and continuous values.
#'
#' This answers a different question from `util_aggregations_allowed()`, and
#' the two must not be confused. This one is a **spell-check**: is the string
#' the caller wrote the name of an aggregation at all, or a typo? The other one
#' is the **semantic** check: is this operator meaningful for *this*
#' `(MARKER, FIGURE)`. A validator needs both, in that order — a helpful error
#' for `"MEDAIN"` ("not a name we know") is not the same as the error for
#' `MODE_LOW` on M-values ("a real name, wrong scale").
#'
#' It exists as its own function because the callers used to obtain the union by
#' asking `util_aggregations_allowed("SIGNAL", "BETA", ...)` — the one pair that
#' happens to admit everything. That worked, but it read as an assumption about
#' the scale of the session, which it never was: a reader on an `MVALUE` run had
#' every reason to think it was a bug.
#'
#' @return character vector of aggregation names.
#' @keywords internal
#' @noRd
util_aggregation_vocabulary <- function() {
  # The union over both classes. SIGNAL/BETA is what widens the set — it is the
  # only combination carrying MODE_LOW/MODE_HIGH — but here it names the axis,
  # not the run: the vocabulary of a language does not depend on the sentence.
  unique(unlist(lapply(c(TRUE, FALSE), function(discrete)
    util_aggregations_allowed("SIGNAL", "BETA", discrete = discrete,
                              default = FALSE))))
}
