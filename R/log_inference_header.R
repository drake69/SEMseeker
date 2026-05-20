#' Log a per-inference-detail journal header
#'
#' Extracted from association_analysis() (was inline at lines 124-134).
#' Emits the JOURNAL banner with the prettified inference_detail row
#' rendered as a kable table.
#'
#' @param inference_detail single-row data.frame.
#' @return Invisibly NULL.
#' @keywords internal
log_inference_header <- function(inference_detail) {
  prettified <- t(inference_detail)
  prettified <- knitr::kable(prettified, format = "simple",
    align = "l", digits = 2, row.names = TRUE)
  prettified <- paste(prettified, collapse = "\n")
  log_event("JOURNAL: ##############################################################################################################")
  log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"),
    " \nStarting association analysis for inference detail:\n", prettified)
  invisible(NULL)
}
