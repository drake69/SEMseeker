#' Coerce a value to logical
#'
#' Converts any representation of TRUE/FALSE (logical, character, or numeric) to
#' a proper R logical value. Returns \code{FALSE} for \code{NA}, empty string,
#' or \code{NULL}.
#'
#' @param x A scalar: logical, character (\code{"TRUE"}, \code{"FALSE"},
#'   \code{"true"}, \code{"T"}, \code{"1"}) or numeric (\code{1}, \code{0}).
#'
#' @return A single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @examples
#' SEMseeker:::boolean_check("TRUE")   # TRUE
#' SEMseeker:::boolean_check("false")  # FALSE
#' SEMseeker:::boolean_check(1)        # TRUE
#' SEMseeker:::boolean_check(NA)       # FALSE
#'
boolean_check <- function(x)
{
  if (is.na(x) || x=="" || is.null(x))
    return(FALSE)

  if (!is.logical(x))
    x = (x == "TRUE" || x == "T" || x == "true" || x == "True" || x == "1" || x == 1)


  return(x)
}
