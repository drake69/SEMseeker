#' Column names of the per-sample statistics sibling (internal)
#'
#' AI-223. These names are the contract between the `SCOPE = SAMPLE` artefacts
#' and the models fitted on them: [sem_study_summary_get()] composes a column per
#' artefact, `association_analysis()` looks its features up by name. Both sides
#' MUST go through these helpers — if they built the strings independently the
#' contract would break at the first divergence.
#'
#' AI-255 removed the file they used to live in (`SAMPLE_STATS_RESULT.csv`); the
#' names outlived it because they name a *feature*, not a file.
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
#' @section Depth-free entry (AI-223 slice 2a):
#' The producer is depth-agnostic by design — the sibling is written ONCE with
#' every requested scope, not once per depth. Called with `depth = NULL` (e.g.
#' `io_scope_name(area = "GENE", subarea = "TSS1500")`) the scope is derived
#' from the pair alone: no area means `SAMPLE`, otherwise `<AREA>[_<SUBAREA>]`.
#' The depth-driven form is what the depth=2/3 consumers still use.
#'
#' @keywords internal
#' @noRd
io_scope_name <- function(depth = NULL, area = NULL, subarea = NULL) {

  if (is.null(depth)) {
    if (is.null(area) || !nzchar(as.character(area)))
      return("SAMPLE")
    parts <- as.character(area)
    if (!is.null(subarea) && nzchar(as.character(subarea)))
      parts <- c(parts, as.character(subarea))
    return(core_name_cleaning(paste(parts, collapse = "_")))
  }

  depth <- suppressWarnings(as.integer(depth))
  if (length(depth) != 1L || is.na(depth)) depth <- 1L

  if (depth <= 1L || is.null(area) || !nzchar(as.character(area)))
    return("SAMPLE")

  parts <- as.character(area)
  if (depth >= 3L && !is.null(subarea) && nzchar(as.character(subarea)))
    parts <- c(parts, as.character(subarea))

  core_name_cleaning(paste(parts, collapse = "_"))
}
