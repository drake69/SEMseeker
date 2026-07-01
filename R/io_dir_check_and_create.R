#' Create a directory path, building any missing intermediate directories
#'
#' Splits \code{baseFolder} into its path components, appends \code{subFolders},
#' and creates each directory level that does not yet exist.  Equivalent to
#' \code{dir.create(path, recursive = TRUE)} but also returns the final
#' normalised absolute path.
#'
#' @param baseFolder Character scalar: root directory path (need not exist yet).
#' @param subFolders Character vector: one or more subdirectory names to append
#'   below \code{baseFolder}.  Each element becomes one level of the hierarchy.
#'
#' @return Character scalar: the normalised absolute path of the deepest
#'   directory created (or already existing).
#'
io_dir_check_and_create <- function(baseFolder, subFolders)
{
  subFolders <- as.vector(vapply(subFolders, as.character, character(1)))
  path <- do.call(file.path, c(list(as.character(baseFolder)), as.list(subFolders)))
  if (!dir.exists(path))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}
