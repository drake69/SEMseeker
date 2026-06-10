#' @title Sensitivity / specificity analysis (deprecated alias)
#' @description Alias deprecato di [diagnostic_performance()]. Il nome
#'   `ss_analysis` era ambiguo (poteva indicare SEMseeker, sensitivity /
#'   specificity, oppure single-sample). Il nuovo nome
#'   `diagnostic_performance()` comunica chiaramente cosa fa: calcolo di
#'   sensitivity, specificity e diagnostic score per le SEM mutations /
#'   lesions.
#'
#'   Sostituire le chiamate `ss_analysis(...)` con
#'   `diagnostic_performance(...)` (stessa firma, stessa semantica). Alla
#'   chiamata viene emesso un warning di deprecation.
#'
#' @param ... Tutti gli argomenti vengono inoltrati a
#'   [diagnostic_performance()] inalterati. Vedi la sua documentazione
#'   per la lista completa.
#' @return Invisibly \code{NULL}. Vedi [diagnostic_performance()].
#' @seealso [diagnostic_performance()]
#' @export
#' @examples
#' # Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
#' # runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
#' invisible(NULL)
ss_analysis <- function(...) {
  .Deprecated("diagnostic_performance", package = "SEMseeker",
              msg = paste0("'ss_analysis' is deprecated. ",
                           "Use 'diagnostic_performance' instead."))
  diagnostic_performance(...)
}
