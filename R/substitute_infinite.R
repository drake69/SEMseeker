substitute_infinite <- function(x) {

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " replacing infinite values with max/min values.")

  x <- as.matrix(x)
  is_inf <- is.infinite(x)
  n_inf_values <- sum(is_inf)

  if (n_inf_values == 0)
    return(as.data.frame(x))

  # Vectorized: replace +Inf/-Inf with max/min finite value (preserving sign)
  max_abs_value <- max(abs(x[is.finite(x)]), na.rm = TRUE)
  x[is_inf] <- sign(x[is_inf]) * max_abs_value

  log_event("JOURNAL: Replaced ", n_inf_values, " infinite values with max/min values.")
  as.data.frame(x)
}
