# ---------------------------------------------------------------------------
# Module-private helpers and tables used by init_env()
#
# Refactor 2026-05-31: collapsed 5 repeated patterns (22 inline
# set_env_variable defaults, 7 inline dir_check_and_create output folders,
# Sys.time-stamped log_event boilerplate, manual NULL/empty filter loop,
# pop-arg-with-default boilerplate) into table-driven plus helper calls.
# Behaviour preserved; only the shape changes.
# ---------------------------------------------------------------------------

# Per-option defaults applied by .apply_defaults() — each entry is a list
# with `value` (mandatory) and optional `choices` (validated by
# set_env_variable when present).
.SS_DEFAULTS <- list(
  verbosity              = list(value = 1),
  q_b_param              = list(value = data.frame("DELTAP_B"=4,"DELTARP_B"=4,"DELTAQ_Q"=4,"DELTARQ_Q"=4)),
  DELTAP_B               = list(value = 4),
  DELTARP_B              = list(value = 4),
  DELTAQ_Q               = list(value = 4),
  DELTARQ_Q              = list(value = 4),
  inpute                 = list(value = "none"),
  plot_format            = list(value = "png"),
  plot_resolution        = list(value = "print"),
  plot_resolution_ppi    = list(value = 600),
  alpha                  = list(value = 0.05),
  sex_chromosome_remove  = list(value = FALSE),
  opencl                 = list(value = FALSE),
  bonferroni_threshold   = list(value = 0.05),
  iqrTimes               = list(value = 3),
  sliding_window_size    = list(value = 11),
  tech                   = list(value = ""),
  genome_build           = list(value = "hg19", choices = c("hg19","hg38","mm10","legacy")),
  showprogress           = list(value = FALSE),
  openai_api_key         = list(value = ""),
  multiple_test_adj      = list(value = "q", choices = c("BY","fdr","BH","bonferroni","q"))
)

# Result-folder subdirectories created under `result_folder`.
.SS_FOLDERS <- c(
  result_folderData       = "Data",
  result_folderChart      = "Chart",
  result_folderInference  = "Inference",
  result_folderPathway    = "Pathway",
  result_folderPhenotype  = "Phenotype",
  result_folderEuristic   = "Euristic",
  session_folder          = "Log"
)

.apply_defaults <- function(arguments, defaults) {
  for (key in names(defaults)) {
    d <- defaults[[key]]
    if (is.null(d$choices))
      arguments <- set_env_variable(arguments, key, d$value)
    else
      arguments <- set_env_variable(arguments, key, d$value, d$choices)
  }
  arguments
}

.log_info <- function(...) {
  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y"), ...)
}

.pop_arg <- function(args, name, default) {
  val <- if (!is.null(args[[name]])) args[[name]] else default
  args[[name]] <- NULL
  list(value = val, args = args)
}

#' init ssEnvonment
#'
#' @param result_folder where result of semseeker will be stored
#' @param maxResources percentage of how many available cores will be used
#'   (default 90 percent, rounded to lowest integer)
#' @param ... additional session options, including:
#'   \describe{
#'     \item{\code{parallel_strategy}}{parallelisation strategy for \pkg{future}:
#'       \code{"sequential"} (default), \code{"multisession"}, \code{"multicore"},
#'       \code{"cluster"}}
#'     \item{\code{genome_build}}{reference genome assembly: \code{"hg19"} (default,
#'       matches Illumina array annotation), \code{"hg38"} (GRCh38, typical for
#'       long-read / Nanopore data), \code{"mm10"} (mouse — requires C-05).
#'       Stored in \code{ssEnv$genome_build} and written to session provenance
#'       metadata (C-06).}
#'     \item{\code{tech}}{override technology auto-detection: \code{"K850"},
#'       \code{"K450"}, \code{"K27"}, \code{"WGBS"}, \code{"LONGREAD"}.
#'       Required for long-read data because LONGREAD cannot be distinguished from
#'       WGBS by probe-ID pattern alone.  Example:
#'       \code{init_env(folder, tech = "LONGREAD", genome_build = "hg38")}}
#'   }
#'
#' @return the working ssEnvonment
init_env <- function(result_folder, maxResources = 90, ...)
{
  gc()
  tryCatch(
    { test_it <- list(...) },
    error = function(cond)  {
      log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Function's arguments must be passed explicitily !")
      log_event(cond)
      stop("Function's arguments must be passed explicitily !")
    }
  )

  withr::local_options(list(digits = 22))

  # suppress warnings messages of packages
  PKGs <- c("future","doRNG","doParallel","progressr","data.table","ggplot2","dplyr",
            "readr","readxl","stringr","tidyr","tibble","purrr","ggpubr","ggrepel","ggsci",
            "foreach","VennDiagram")
  tt <- lapply(PKGs, suppressWarnings(suppressMessages))
  tt <- lapply(PKGs, suppressPackageStartupMessages)

  arguments <- list(...)
  if (length(arguments) > 0) {
    # remove empty character items but preserve logical/numeric
    arguments <- lapply(arguments, function(x) if (is.character(x)) gsub(" ", "", x) else x)
    arguments <- lapply(arguments, function(x) if (is.character(x)) x[x != ""] else x)
    arguments <- arguments[sapply(arguments, function(x) length(x) > 0)]
    arguments <- arguments[sapply(arguments, function(x) !is.null(x))]
  }

  arguments[["areas_selection"]] <- NULL

  popped <- .pop_arg(arguments, "start_fresh", FALSE)
  start_fresh <- popped$value
  arguments   <- popped$args

  if (start_fresh) {
    unlink(result_folder, recursive = TRUE, force = TRUE)
    ssEnv <- list()
  } else if (dir.exists(result_folder)) {
    ssEnv <- get_session_info(result_folder)
  } else {
    ssEnv <- list()
  }

  if (is.null(ssEnv$session_id))
    ssEnv$session_id <- 0
  else
    ssEnv$session_id <- ssEnv$session_id + 1
  ssEnv$session_folder <- dir_check_and_create(result_folder, c("Log"))
  update_session_info(ssEnv)

  ssEnv$seed <- 7658776

  arguments <- .apply_defaults(arguments, .SS_DEFAULTS)

  if (!is.null(ssEnv$openai_api_key) && nzchar(ssEnv$openai_api_key))
    message("SEMseeker: set OPENAI_API_KEY in your environment to enable OpenAI features.")

  # tech × genome_build validation is centralised in semseeker() (the public
  # dispatcher); init_env() no longer duplicates the check.

  original_colors <- c('#b9e192', '#b3c7f7', '#f8b8d0','#f194b8', '#ffefb6', '#cfebb6','#b9ef92')
  original_colors <- rep(original_colors, 2)
  arguments <- set_env_variable(arguments, "color_palette", original_colors)
  darker_colors <- grDevices::adjustcolor(original_colors, alpha.f = 0.5)
  darker_colors <- c("blue","red","purple","green","yellow","orange","brown")
  arguments <- set_env_variable(arguments, "color_palette_darker", darker_colors)
  arguments <- set_env_variable(arguments, "cluster_workers", NULL)

  model_metrics <- toupper(as.vector(SEMseeker::metrics_properties$Metric))
  arguments <- set_env_variable(arguments, "model_metrics", model_metrics)

  ssEnv <- get_session_info()

  popped  <- .pop_arg(arguments, "dry_run", FALSE)
  dry_run <- popped$value
  arguments <- popped$args
  if (dry_run)
    ssEnv$verbosity <- 4

  tmp <- tempdir()
  .log_info(" data will saved in this folder:", result_folder)
  ssEnv$temp_folder   <- paste(tmp, "/semseeker/",
                               stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),
                               sep = "")
  ssEnv$result_folder <- result_folder

  for (key in names(.SS_FOLDERS))
    ssEnv[[key]] <- dir_check_and_create(result_folder, .SS_FOLDERS[[key]])

  random_file_name <- paste(stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"), ".log", sep = "")

  # Skip sink when running inside a callr child process — callr already
  # redirects stdout/stderr to log.txt; opening a second sink fights with
  # the redirect and can cause silent crashes on large I/O.
  if (!identical(Sys.getenv("SEMSEEKER_CHILD"), "1")) {
    if (sink.number() != 0)
      sink(NULL)
    file_name <- paste(as.character(Sys.info()["nodename"]), "_session_output.log", sep = "")
    sink(file.path(ssEnv$session_folder, file_name), split = TRUE, append = TRUE)
  }

  foreachIndex <- 0

  arguments <- keys_create(ssEnv, arguments)
  ssEnv <- get_session_info()

  ssEnv$functionToExport <- c("analyze_single_sample","deltar_single_sample",
                              "dump_sample_as_bed_file", "delta_single_sample",
                              "dir_check_and_create", "file_path_build",
                              "analyze_single_sample_both", "sort_by_chr_and_start",
                              "test_match_order", "lesions_get", "mutations_get")

  if (ssEnv$showprogress) {
    handler_settings <- progressr::handlers()
    if (!(exists("cli", mode = "function", inherits = TRUE))) {
      if (!testthat::is_testing())
        if (!("cli" %in% handler_settings$handler)) {
          progressr::handlers(global = TRUE)
          progressr::handlers("cli")
        }
    }
  }

  arguments <- set_env_variable(arguments, "maxResources", maxResources)
  arguments <- set_env_variable(arguments, "parallel_strategy", "sequential")
  parallel_session()
  ssEnv <- get_session_info()

  .log_info(" I will focus on:",
            paste(unique(ssEnv$keys_markers_figures$MARKER), collapse = " ", sep = " "),
            " due to ",
            paste(unique(ssEnv$keys_markers_figures$FIGURE), collapse = " ", sep = " "),
            " of ",
            paste(unique(ssEnv$keys_areas_subareas_markers_figures$AREA), collapse = " ", sep = " "))

  # drop NULL / empty character() arguments using Filter (replaces an
  # explicit manual index loop)
  arguments <- Filter(function(x) !is.null(x) && !identical(x, character(0)), arguments)

  if (length(arguments) != 0) {
    .log_info(" This options are not recognized: ",
              paste(arguments, collapse = " ", sep = " "))
    stop("ERROR: This options are not recognized: ",
         paste(arguments, collapse = " ", sep = " "))
  }

  if (dry_run) {
    knitr::kable(as.data.frame(ssEnv$keys_areas_subareas_markers_figures),
                 format = "pipe", caption = "Selection:")
    message(ssEnv$keys_areas_subareas_markers_figures)
    stop("INFO: Dry run is requested. Exiting now.")
  }

  update_session_info(ssEnv)
  return(ssEnv)
}
