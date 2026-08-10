#' Column names of the per-sample statistics sibling (internal)
#'
#' AI-223. The sibling artefact (`SAMPLE_STATS_RESULT.csv`) is the contract
#' between production (annotation + aggregation, inside `semseeker()`) and
#' consumption (`association_analysis()`). Both sides MUST compose its column
#' names through these helpers: if producer and consumer built the strings
#' independently, the contract would break at the first divergence.
#'
#' Naming, all uppercase via [core_name_cleaning()]:
#' \itemize{
#'   \item scope — `SAMPLE` (whole sample, no genomic filter), `<AREA>`, or
#'     `<AREA>_<SUBAREA>`;
#'   \item signal descriptors (marker-independent) — `<SCOPE>_<STAT>`, e.g.
#'     `SAMPLE_MEDIAN`, `GENE_BODY_IQR`;
#'   \item burden per marker — `<SCOPE>_<MARKER>_<FIGURE>`, e.g.
#'     `SAMPLE_MUTATIONS_HYPER`.
#' }
#'
#' @keywords internal
#' @noRd
io_scope_name <- function(depth = 1L, area = NULL, subarea = NULL) {
  depth <- suppressWarnings(as.integer(depth))
  if (length(depth) != 1L || is.na(depth)) depth <- 1L

  if (depth <= 1L || is.null(area) || !nzchar(as.character(area)))
    return("SAMPLE")

  parts <- as.character(area)
  if (depth >= 3L && !is.null(subarea) && nzchar(as.character(subarea)))
    parts <- c(parts, as.character(subarea))

  core_name_cleaning(paste(parts, collapse = "_"))
}

#' @keywords internal
#' @noRd
io_stat_colname <- function(scope, stat) {
  if (missing(scope) || is.null(scope) || !nzchar(scope))
    stop("io_stat_colname(): 'scope' is required (use io_scope_name()).")
  core_name_cleaning(paste(scope, stat, sep = "_"))
}

#' @keywords internal
#' @noRd
io_burden_colname <- function(scope, marker, figure) {
  if (missing(scope) || is.null(scope) || !nzchar(scope))
    stop("io_burden_colname(): 'scope' is required (use io_scope_name()).")
  core_name_cleaning(paste(scope, marker, figure, sep = "_"))
}

#' Statistics emitted for a signal scope
#'
#' MODE_LOW / MODE_HIGH only make sense on a bounded, bimodal scale (beta in
#' [0,1], the two peaks around 0 and 1). On M-values the distribution is
#' roughly gaussian and the two-mode split carries no meaning, so the columns
#' are omitted rather than filled with a meaningless number.
#'
#' @keywords internal
#' @noRd
io_signal_stats <- function(beta = TRUE) {
  base <- c("MEDIAN", "MEAN", "VARIANCE", "IQR", "N_PROBES")
  if (isTRUE(beta)) c(base, "MODE_LOW", "MODE_HIGH") else base
}
