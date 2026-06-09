#' @title Association Cross Studies Overlaps (deprecated alias)
#' @description Alias deprecato di [inter_study_association_overlaps()]. Naming uniformato a
#'   `<scope>_<object>_<operation>` (es. `intra_study_*`,
#'   `inter_study_*`). Sostituire le chiamate `association_cross_studies_overlaps(...)` con
#'   `inter_study_association_overlaps(...)` (stessa firma). Alla chiamata viene emesso un warning
#'   di deprecation.
#' @param ... Argomenti inoltrati a `inter_study_association_overlaps()`.
#' @return Invisibly `NULL`. Vedi [inter_study_association_overlaps()].
#' @seealso [inter_study_association_overlaps()]
#' @export
association_cross_studies_overlaps <- function(...) {
  .Deprecated("inter_study_association_overlaps", package = "SEMseeker",
              msg = paste0("'association_cross_studies_overlaps' is deprecated. ",
                           "Use 'inter_study_association_overlaps' instead."))
  inter_study_association_overlaps(...)
}
