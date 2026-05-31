test_that("log-event", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90,
    figures = "HYPER", markers = "DELTAS", areas = "GENE")

  file_name <- paste(as.character(Sys.info()["nodename"]),"_session_output.log", sep="")

  # log event
  SEMseeker:::log_event("GEOquery package is not installed. Please install pathfindR package to use this function")

  # get session_log file number of row
  session_log <- read.table(file.path(ssEnv$session_folder, file_name), sep = "\t", header = TRUE, quote = "")

  # check if the log event is in the session_log file
  expect_true(any(grepl("GEOquery package is not installed. Please install pathfindR package to use this function", session_log[1])))

  # add event
  SEMseeker:::log_event("Phenolizer: phenotype.hpoa file not found")

  # get session_log file number of row
  session_log <- read.table(file.path(ssEnv$session_folder, file_name), sep = "\t", header = TRUE, quote = "")

  # check row are increased
  testthat::expect_true(nrow(session_log) > 15)

  ####################################################################################
  SEMseeker:::close_env()

  session_log_closed <- read.table(file.path(ssEnv$session_folder, file_name), sep = "\t", header = TRUE, quote = "")
  testthat::expect_true(nrow(session_log) == nrow(session_log_closed) - 2)

  # check after closing session log is preserved when not starting fresh
  ssEnv <- SEMseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "HYPER", markers = "DELTAS", areas = "GENE", start_fresh = FALSE)
  session_log <- read.table(file.path(ssEnv$session_folder, file_name), sep = "\t", header = TRUE, quote = "")
  testthat::expect_true(nrow(session_log) > nrow(session_log_closed))

  unlink(tempFolder, recursive = TRUE)

})
