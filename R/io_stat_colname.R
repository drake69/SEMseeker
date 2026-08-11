#' Compose a signal-descriptor column name of the sibling artefact
#'
#' `<SCOPE>_<STAT>`, e.g. `SAMPLE_MEDIAN` or `GENE_BODY_IQR`. See
#' [io_scope_name()] for the naming contract this belongs to.
#'
#' @param scope scope name, as returned by [io_scope_name()].
#' @param stat one of the statistics declared by [io_signal_stats()].
#' @keywords internal
#' @noRd
io_stat_colname <- function(scope, stat) {
  if (missing(scope) || is.null(scope) || !nzchar(scope))
    stop("io_stat_colname(): 'scope' is required (use io_scope_name()).")
  core_name_cleaning(paste(scope, stat, sep = "_"))
}
