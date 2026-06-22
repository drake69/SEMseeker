#' @title Pathway Cross Subsamples Overlaps (deprecated alias)
#' @description Alias deprecato di [intra_study_pathway_subsamples_overlaps()]. Naming uniformato a
#'   `<scope>_<object>_<operation>` (es. `intra_study_*`,
#'   `inter_study_*`). Sostituire le chiamate `pathway_cross_subsamples_overlaps(...)` con
#'   `intra_study_pathway_subsamples_overlaps(...)` (stessa firma). Alla chiamata viene emesso un warning
#'   di deprecation.
#' @param ... Argomenti inoltrati a `intra_study_pathway_subsamples_overlaps()`.
#' @return Invisibly `NULL`. Vedi [intra_study_pathway_subsamples_overlaps()].
#' @seealso [intra_study_pathway_subsamples_overlaps()]
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
#' invisible(NULL)
pathway_cross_subsamples_overlaps <- function(...) {
  .Deprecated("intra_study_pathway_subsamples_overlaps", package = "SEMseeker",
              msg = paste0("'pathway_cross_subsamples_overlaps' is deprecated. ",
                           "Use 'intra_study_pathway_subsamples_overlaps' instead."))
  intra_study_pathway_subsamples_overlaps(...)
}
