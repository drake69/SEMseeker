#' @title Association Cross Studies Meta Analysis (deprecated alias)
#' @description Alias deprecato di [inter_study_association_meta_analysis()]. Naming uniformato a
#'   `<scope>_<object>_<operation>` (es. `intra_study_*`,
#'   `inter_study_*`). Sostituire le chiamate `association_cross_studies_meta_analysis(...)` con
#'   `inter_study_association_meta_analysis(...)` (stessa firma). Alla chiamata viene emesso un warning
#'   di deprecation.
#' @param ... Argomenti inoltrati a `inter_study_association_meta_analysis()`.
#' @return Invisibly `NULL`. Vedi [inter_study_association_meta_analysis()].
#' @seealso [inter_study_association_meta_analysis()]
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
#' invisible(NULL)
association_cross_studies_meta_analysis <- function(...) {
  .Deprecated("inter_study_association_meta_analysis", package = "SEMseeker",
              msg = paste0("'association_cross_studies_meta_analysis' is deprecated. ",
                           "Use 'inter_study_association_meta_analysis' instead."))
  inter_study_association_meta_analysis(...)
}
