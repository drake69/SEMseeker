core_close_env <- function()
{

  ssEnv <- core_get_session_info()

  if (ssEnv$showprogress)
    progressr::handlers()

  # build all the folder tree in result_folder

  core_remove_empty_folders(ssEnv$result_folder)

  future::plan( future::sequential)
  unlink(ssEnv$temp_folder,recursive=TRUE)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Job Completed !")
  core_log_event("DEBUG: --------------------------------------------------------------")

}


core_remove_empty_folders <- function(path) {
  files <- list.files(path, full.names = TRUE)

  for (file in files) {
    if (file.info(file)$isdir) {
      core_remove_empty_folders(file)  # Recursively check subdirectories

      # After the check, remove the directory if empty
      if (length(list.files(file)) == 0) {
        unlink(file, recursive = TRUE)
        # cat("Removed empty directory:", file, "\n")
      }
    }
  }
}
