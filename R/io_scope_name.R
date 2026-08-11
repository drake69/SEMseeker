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
#'     `SAMPLE_MEDIAN`, `GENE_BODY_IQR`, composed by [io_stat_colname()];
#'   \item burden per marker — `<SCOPE>_<MARKER>_<FIGURE>`, e.g.
#'     `SAMPLE_MUTATIONS_HYPER`, composed by [io_burden_colname()].
#' }
#'
#' The set of statistics a scope carries is declared by [io_signal_stats()] and
#' computed by [util_signal_descriptors()].
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
