#' Compose a per-marker burden column name of the sibling artefact
#'
#' `<SCOPE>_<MARKER>_<FIGURE>`, e.g. `SAMPLE_MUTATIONS_HYPER`. See
#' [io_scope_name()] for the naming contract this belongs to.
#'
#' @param scope scope name, as returned by [io_scope_name()].
#' @param marker marker name.
#' @param figure figure name (`HYPER` / `HYPO`).
#' @keywords internal
#' @noRd
io_burden_colname <- function(scope, marker, figure) {
  if (missing(scope) || is.null(scope) || !nzchar(scope))
    stop("io_burden_colname(): 'scope' is required (use io_scope_name()).")
  core_name_cleaning(paste(scope, marker, figure, sep = "_"))
}
