recover <- function(result_folder, maxResources = 90,
  parallel_strategy  = "multicore",start_fresh = FALSE, ...)
{

  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will try to recover the analysys for project \n in ", ssEnv$result_folderData)

  study_summary <-   study_summary_get()

  create_position_pivots(study_summary,ssEnv$keys_markers_figures)
  deltaX_get()
  study_summary_total()
  annotate_position_pivots()
}
