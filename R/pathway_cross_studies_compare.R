#' @title Pathway Cross Studies Compare (deprecated alias)
#' @description Alias deprecato di [inter_study_pathway_compare()]. Naming uniformato a
#'   `<scope>_<object>_<operation>` (es. `intra_study_*`,
#'   `inter_study_*`). Sostituire le chiamate `pathway_cross_studies_compare(...)` con
#'   `inter_study_pathway_compare(...)` (stessa firma). Alla chiamata viene emesso un warning
#'   di deprecation.
#' @param ... Argomenti inoltrati a `inter_study_pathway_compare()`.
#' @return Invisibly `NULL`. Vedi [inter_study_pathway_compare()].
#' @seealso [inter_study_pathway_compare()]
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
#' invisible(NULL)
pathway_cross_studies_compare <- function(...) {
  .Deprecated("inter_study_pathway_compare", package = "SEMseeker",
              msg = paste0("'pathway_cross_studies_compare' is deprecated. ",
                           "Use 'inter_study_pathway_compare' instead."))
  inter_study_pathway_compare(...)
}
