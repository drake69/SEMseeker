#' Convert M-values to beta values
#'
#' Applies the standard transformation \eqn{\beta = 2^M / (1 + 2^M)} element-wise.
#' Handles numeric matrices and numeric columns of data frames. Non-numeric
#' columns (e.g. annotation columns in coordinate-format data frames) are
#' preserved unchanged.
#'
#' @param x A numeric matrix, numeric vector, or data frame containing M-values.
#' @param coord_cols Optional character vector of column names to preserve
#'   unchanged (e.g. \code{c("CHR","START","END")}). If \code{NULL} (default),
#'   all numeric columns are converted; all non-numeric columns are preserved.
#'
#' @return An object of the same class as \code{x} with M-values converted to
#'   beta values.
#' @keywords internal
mvalue_to_beta <- function(x, coord_cols = NULL) {
  m_to_b <- function(m) 2^m / (1 + 2^m)

  if (is.matrix(x) || is.numeric(x)) {
    return(m_to_b(x))
  }

  if (is.data.frame(x)) {
    preserve <- if (is.null(coord_cols)) {
      vapply(x, function(col) !is.numeric(col), logical(1L))
    } else {
      colnames(x) %in% coord_cols
    }
    for (j in which(!preserve)) {
      x[[j]] <- m_to_b(x[[j]])
    }
    return(x)
  }

  stop("mvalue_to_beta(): unsupported input class '", class(x)[1L], "'.")
}
