#' @title Gene impact analysis (deprecated, use enrichment_analysis)
#' @description Alias deprecato di [enrichment_analysis()]. La funzione era
#'   originariamente esposta come `gene_impact_analysis` ma il nome non
#'   riflette accuratamente cosa fa (= dispatcher di pathway/phenotype
#'   enrichment). Mantenuta come alias per backward-compatibility col codice
#'   utente esistente; emette un warning di deprecation alla chiamata.
#'
#'   Sostituire le chiamate `gene_impact_analysis(...)` con
#'   `enrichment_analysis(...)` (stessa firma, stessa semantica).
#'
#' @param ... Argomenti forwardati a `enrichment_analysis()`.
#' @return Invisibly NULL. Vedi `enrichment_analysis()`.
#' @seealso [enrichment_analysis()]
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
#' invisible(NULL)
gene_impact_analysis <- function(...) {
  .Deprecated("enrichment_analysis", package = "SEMseeker",
              msg = paste0("'gene_impact_analysis' is deprecated. ",
                           "Use 'enrichment_analysis' instead."))
  enrichment_analysis(...)
}
