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
dir_check_and_create <- function  (baseFolder, subFolders)
{


  folderSep <- as.character(.Platform$file.sep)
  parts <- unlist(strsplit(as.character(baseFolder), folderSep))
  # do.call(file.path, as.list(c(parts[seq_along(parts) - 1], subFolders)))
  subFolders <- as.vector(sapply(subFolders, as.character))
  subFolders <-c(parts[seq_along(parts)], subFolders)

  for( y in seq_along(subFolders))
  {

    subFolder <- subFolders[y]

    # if(subFolder=="c(\"Reference\", \"Control\", \"Case\")")
    #   #

    # # browser
    if(!exists("subF"))
      subF <- file.path(subFolder)
    else
    {
      subF <- file.path(subF, subFolder)
      if(!dir.exists(subF))
        dir.create(subF, recursive = FALSE)
    }
  }


  res <- file.path(as.character(subF))

  # transform to absolute path
  res <- file.path(normalizePath(res))

  return(res)
}
