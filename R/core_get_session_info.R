core_get_session_info <- function(result_folder=NULL)
{

  # try from environment
  ssEnv <- .pkgglobalenv$ssEnv

  session_folder <- file.path(result_folder,"/Log/session_info.rds")
  # try from file
  if (is.null(ssEnv) | length(ssEnv)==0)
    if((!is.null(result_folder) & file.exists(session_folder)))
      ssEnv <- readRDS( file.path(session_folder))

  # try from doFuture
  # if (is.null(ssEnv) | length(ssEnv)==0)
  #   ssEnv <- getOption("ssEnv")

  if (is.null(ssEnv) | length(ssEnv)==0) {
    if (!is.null(result_folder))
      return(list())
    stop("ERROR: core_get_session_info called without result folder!")
  }

  return(ssEnv)

}
