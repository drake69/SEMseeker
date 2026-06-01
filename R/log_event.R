log_event <- function(...)
{
  # Robustness: log_event must be callable even when no session has been
  # initialised (e.g. area_granges_build() invoked standalone, outside of a
  # full semseeker() run). Silently no-op instead of crashing on missing ssEnv.
  ssEnv <- tryCatch(get_session_info(), error = function(e) NULL)
  if (is.null(ssEnv) || length(ssEnv) == 0L || is.null(ssEnv$session_folder))
    return(invisible(NULL))

  # append log_event to log file
  log_events <- list(...)
  log_event_to_save <- ""
  for (i in seq_along(log_events))
  {
    log_event_to_save <- paste0(log_event_to_save, log_events[i], sep =" ")
  }

  log_event_to_save <- gsub("  "," ", log_event_to_save)
  log_event_to_save <- gsub("  "," ", log_event_to_save)
  log_event_to_save <- gsub("  "," ", log_event_to_save)
  log_event_to_save <- gsub("  "," ", log_event_to_save)

  # Auto-augment DEBUG messages with memory + CPU usage. Calling `gc()` and
  # `ps` is cheap enough for debug-level granularity but is intentionally NOT
  # done for INFO/WARNING/ERROR to keep hot paths fast and logs compact.
  if (grepl("^DEBUG", log_event_to_save)) {
    mem_mb <- tryCatch(
      sum(gc(verbose = FALSE, full = FALSE)[, "(Mb)"]),
      error = function(e) NA_real_
    )
    cpu_pct <- tryCatch(
      as.numeric(trimws(system(
        sprintf("ps -p %d -o %%cpu=", Sys.getpid()),
        intern = TRUE)[1])),
      error = function(e) NA_real_
    )
    log_event_to_save <- paste0(log_event_to_save,
      " [mem=", round(mem_mb, 1), "MB cpu=", round(cpu_pct, 1), "%]")
  }

  file_name <- paste(as.character(Sys.info()["nodename"]),"_session_output.log", sep="")
  log_file <- file.path(ssEnv$session_folder,file_name)

  log_event_to_print <- log_event_to_save
  if (grepl("^BANNER", log_event_to_save))
  {
    log_event_to_print_dash <- "##############################################################################"
    log_event_to_print <- paste0( log_event_to_print_dash, "\n",
      gsub("BANNER: ", "", log_event_to_print), "\n",
      log_event_to_print_dash, '\n')
  }

  if (grepl("ERROR:", log_event_to_save))
  {
    log_event_to_print_dash <- "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    log_event_to_print <- paste0( log_event_to_print_dash, "\n",  log_event_to_print  , "\n",log_event_to_print_dash)
  }

  if(is.null(ssEnv$verbosity))
    verbosity <- 4
  else
    verbosity <- as.numeric(ssEnv$verbosity)

  cat(log_event_to_save, "\n", file = log_file, append = TRUE)

  file_name <- "console_session_output.log"
  log_file <- file.path(ssEnv$session_folder,file_name)
  cat(log_event_to_save, "\n", file = log_file, append = TRUE)

  if(grepl("^JOURNAL", log_event_to_save))
  {
    file_name <- "lab_journal.log"
    journal_file <- file.path(ssEnv$session_folder,file_name)
    log_event_to_save <- gsub("^JOURNAL: ", "", log_event_to_save)
    cat(log_event_to_save, "\n", file = journal_file, append = TRUE)
  }

  # if verbosity is 1, only print ERROR messages
  if (verbosity == 1 && grepl("^ERROR", log_event_to_save) || grepl("^BANNER", log_event_to_save) && !testthat::is_testing())
  {
    message(log_event_to_print)
    log_event_to_print <- ""
  }

  # if verbosity is 2, print ERROR and WARNING messages
  if (verbosity == 2 && (grepl("^ERROR", log_event_to_save) || grepl("^WARNING", log_event_to_save)|| grepl("^BANNER", log_event_to_save)) && !testthat::is_testing())
  {
    message(log_event_to_print)
    log_event_to_print <- ""
  }

  # if verbosity is 3, print ERROR, WARNING and INFO messages
  if (verbosity == 3 && (grepl("^ERROR", log_event_to_save) || grepl("^WARNING", log_event_to_save) || grepl("^INFO", log_event_to_save) || grepl("^BANNER", log_event_to_save)) && !testthat::is_testing())
  {
    message(log_event_to_print)
    log_event_to_print <- ""
  }

  # if verbosity is 4, print ERROR, WARNING, INFO and DEBUG messages
  if (verbosity == 4 && !testthat::is_testing())
  {
    message(log_event_to_print)
    log_event_to_print <- ""
  }

}
