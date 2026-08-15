#' Column names of the per-sample features (internal)
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
#' The name of a scope is the region class it covers — `GENE_TSS1500` — or
#' `SAMPLE` when there is no restriction at all. It used to be derived from
#' `depth`, with `depth <= 1` meaning SAMPLE, `2` meaning the area alone and `3`
#' the area with its subarea. That mapping is gone with depth itself: the class
#' is said by the pair, and a number that ordered pairs on a line could only lie
#' about the ones that are not comparable.
#'
#' @param area,subarea the region class; both `NULL` gives `"SAMPLE"`.
#' @return the scope name, uppercase via [core_name_cleaning()].
#' @keywords internal
#' @noRd
io_scope_name <- function(area = NULL, subarea = NULL) {

  if (is.null(area) || !nzchar(as.character(area)))
    return("SAMPLE")

  parts <- as.character(area)
  if (!is.null(subarea) && nzchar(as.character(subarea)))
    parts <- c(parts, as.character(subarea))

  core_name_cleaning(paste(parts, collapse = "_"))
}
