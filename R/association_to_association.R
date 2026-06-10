#' @title Association To Association (deprecated alias)
#' @description Alias deprecato di [intra_study_association_replication()]. Naming uniformato a
#'   `<scope>_<object>_<operation>` (es. `intra_study_*`,
#'   `inter_study_*`). Sostituire le chiamate `association_to_association(...)` con
#'   `intra_study_association_replication(...)` (stessa firma). Alla chiamata viene emesso un warning
#'   di deprecation.
#' @param ... Argomenti inoltrati a `intra_study_association_replication()`.
#' @return Invisibly `NULL`. Vedi [intra_study_association_replication()].
#' @seealso [intra_study_association_replication()]
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
#' invisible(NULL)
association_to_association <- function(...) {
  .Deprecated("intra_study_association_replication", package = "SEMseeker",
              msg = paste0("'association_to_association' is deprecated. ",
                           "Use 'intra_study_association_replication' instead."))
  intra_study_association_replication(...)
}
