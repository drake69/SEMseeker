#' Find gene sets unique to each group
#'
#' Given a named list where each element is a character vector of gene set
#' identifiers, returns only those identifiers that appear exclusively in one
#' group and not in any other.
#'
#' @param split_list A named list of character vectors (one per group/category).
#'
#' @return A named list of the same structure, containing only the gene sets
#'   that are unique to each group.  Groups with no unique sets are dropped.
#'
find_unique_gene_sets <- function(split_list) {
  unique_sets <- list()
  keys <- names(split_list)

  for (k in seq_along(keys)) {
    key <- keys[k]
    current_sets <- split_list[[key]]
    other_keys <- setdiff(keys, key)
    other_sets <- unlist(split_list[other_keys], use.names = FALSE)
    unique_to_current <- setdiff(current_sets, other_sets)
    unique_sets[[key]] <- unique_to_current
  }

  # remove empty sets
  unique_sets <- unique_sets[sapply(unique_sets, length) > 0]
  return(unique_sets)
}

#' Wrap long strings to a fixed width
#'
#' Wraps each string in \code{x} to at most \code{len} characters per line,
#' collapsing the result with newlines.  Useful for plot axis labels.
#'
#' @param x Character vector of strings to wrap.
#' @param len Integer maximum line width in characters.
#'
#' @return Character vector of the same length as \code{x}, with long strings
#'   broken into multiple lines separated by \code{\\n}.
#'
wrap_it <- function(x, len)
{
  sapply(x, function(y) paste0(strwrap(y, len),
    collapse = "\n"),
    USE.NAMES = FALSE)
}
