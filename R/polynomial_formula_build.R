#' Build a polynomial regression formula with covariate interactions
#'
#' Pure helper: no I/O, no side effects. Constructs a formula of the form
#' \code{y ~ I(x^1) + I(x^2) + ... + I(x^degree) + I(x^1):cov1 + I(x^1):cov2 + ...}
#'
#' @param dependent_variable  Name of the response column (string).
#' @param independent_variable Name of the predictor column (string).
#' @param degree              Polynomial degree (positive integer).
#' @param covariates          Character vector of covariate column names (may be empty).
#' @return A \code{formula} object.
#'
polynomial_formula_build <- function(dependent_variable, independent_variable, degree, covariates) {
  polynomial_terms <- paste0("I(", independent_variable, "^", seq_len(degree), ")")
  x_part <- paste(polynomial_terms, collapse = " + ")

  if (length(covariates) > 0) {
    interaction_terms <- unlist(lapply(polynomial_terms, function(pt) paste0(pt, ":", covariates)))
    interaction_part  <- paste(interaction_terms, collapse = " + ")
    formula_string    <- paste(dependent_variable, "~", x_part, "+", interaction_part)
  } else {
    formula_string <- paste(dependent_variable, "~", x_part)
  }

  stats::as.formula(formula_string)
}
