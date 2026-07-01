#' Title
#'
#' @param x vector to compare
#' @param y vector to compare
#'
#' @return true if the order matches otherwise is false

#'
test_match_order <- function(x, y) {


  if(is.null(x) | is.null(y))
    return(FALSE)

  if (all(x == y)) {
    return(TRUE)
  }
  # log_event('Perfect match in same order')

  if (!all(x == y) && all(sort(x) == sort(y))) {
    return(FALSE)
  }
  # log_event('Perfect match in wrong order')

  if (!all(x == y) && !all(sort(x) == sort(y))) {
    return(FALSE)
  }
  # log_event('No match')
}


test_match_order_by_rownames <- function(x, y) {


  if(is.null(x) | is.null(y))
    return(FALSE)

  if (all(rownames(x) == rownames(y))) {
    return(TRUE)
  }
  return(FALSE)
  # log_event('No match')
}
