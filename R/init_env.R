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
    {
      test_it <- list(...)
    },
    error = function(cond)  {
      log_event ("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Function's arguments must be passed explicitily !")
      log_event(cond)
      stop("Function's arguments must be passed explicitily !")
    }
  )

  # set digits to 22
  withr::local_options(list(digits = 22))

  # suppress warnings messages of packages
  PKGs<- c("future","doRNG","doParallel","progressr","data.table","ggplot2","dplyr",
    "readr","readxl","stringr","tidyr","tibble","purrr","ggpubr","ggrepel","ggsci","foreach","VennDiagram")
  tt <- lapply(PKGs, suppressWarnings(suppressMessages))
  tt <- lapply(PKGs, suppressPackageStartupMessages)


  arguments <- list(...)
  # check if optional arguments are passed
  if(length(arguments) == 0)
  {
    arguments <- list()
  }
  else
  {
    # remove all empty items from arguments (only apply gsub to character, preserve logical/numeric types)
    arguments <- lapply(arguments, function(x) if(is.character(x)) gsub(" ", "", x) else x)
    arguments <- lapply(arguments, function(x) if(is.character(x)) x[x!=""] else x)
    arguments <- arguments[sapply(arguments, function(x) length(x) > 0)]
    arguments <- arguments[sapply(arguments, function(x) !is.null(x))]
    # arguments <- arguments[sapply(arguments, function(x) !is.na(x))]
  }

  arguments[["areas_selection"]] <- NULL


  start_fresh <- FALSE
  if(!is.null(arguments[["start_fresh"]]))
    start_fresh <- arguments$start_fresh
  arguments[["start_fresh"]] <- NULL

  if(start_fresh)
  {
    unlink(result_folder, recursive = TRUE, force = TRUE)
    ssEnv <- list()
  }
  else if(dir.exists(result_folder))
    ssEnv <- get_session_info(result_folder)
  else
    ssEnv <- list()

  if(is.null(ssEnv$session_id))
    ssEnv$session_id <- 0
  else
    ssEnv$session_id <- ssEnv$session_id + 1
  ssEnv$session_folder <-  dir_check_and_create(result_folder,c("Log"))
  update_session_info(ssEnv)

  ssEnv$seed <- 7658776


  arguments <- set_env_variable(arguments,"verbosity",1)
  arguments <- set_env_variable(arguments,"q_b_param",data.frame("DELTAP_B"=4,"DELTARP_B"=4,"DELTAQ_Q"=4,"DELTARQ_Q"=4))
  arguments <- set_env_variable(arguments,"DELTAP_B",4)
  arguments <- set_env_variable(arguments,"DELTARP_B",4)
  arguments <- set_env_variable(arguments,"DELTAQ_Q",4)
  arguments <- set_env_variable(arguments,"DELTARQ_Q",4)

  arguments <- set_env_variable(arguments,"inpute","none")
  arguments <- set_env_variable(arguments,"plot_format","png")
  arguments <- set_env_variable(arguments,"plot_resolution","print")
  arguments <- set_env_variable(arguments,"plot_resolution_ppi",600)
  arguments <- set_env_variable(arguments,"alpha",0.05)
  arguments <- set_env_variable(arguments,"sex_chromosome_remove",FALSE)
  arguments <- set_env_variable(arguments,"opencl",FALSE)
  arguments <- set_env_variable(arguments,"bonferroni_threshold",0.05)
  arguments <- set_env_variable(arguments,"iqrTimes",3)
  arguments <- set_env_variable(arguments,"sliding_window_size",11)
  arguments <- set_env_variable(arguments,"tech","")
  arguments <- set_env_variable(arguments,"genome_build","hg19",
    c("hg19","hg38","mm10","legacy"))
  arguments <- set_env_variable(arguments,"showprogress",FALSE)
  arguments <- set_env_variable(arguments,"openai_api_key","")
  arguments <- set_env_variable(arguments,"multiple_test_adj","q", c("BY", "fdr","BH","bonferroni","q"))

  if (!is.null(ssEnv$openai_api_key) && nzchar(ssEnv$openai_api_key))
    message("SEMseeker: set OPENAI_API_KEY in your environment to enable OpenAI features.")

  # tech × genome_build validation is centralised in semseeker() (the public
  # dispatcher); init_env() no longer duplicates the check.

  original_colors <- c('#b9e192', '#b3c7f7', '#f8b8d0','#f194b8', '#ffefb6', '#cfebb6','#b9ef92')
  original_colors <- rep(original_colors, 2)
  # original_colors <- khroma::color("bright", n = 20)
  arguments <- set_env_variable(arguments,"color_palette",original_colors)
  darker_colors <- grDevices::adjustcolor(original_colors, alpha.f = 0.5)
  darker_colors <- c("blue","red","purple","green","yellow","orange","brown")
  arguments <- set_env_variable(arguments,"color_palette_darker",darker_colors)
  arguments <- set_env_variable(arguments,"cluster_workers",NULL)

  model_metrics <- toupper(as.vector(SEMseeker::metrics_properties$Metric))
  arguments <- set_env_variable(arguments,"model_metrics",model_metrics)

  # get ssEnv
  ssEnv <- get_session_info()

  dry_run <- FALSE
  if(!is.null(arguments[["dry_run"]]))
    dry_run <- arguments$dry_run
  arguments[["dry_run"]] <- NULL
  if(dry_run)
    ssEnv$verbosity <- 4

  tmp <- tempdir()
  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y")," data will saved in this folder:", result_folder)
  ssEnv$temp_folder <-  paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  ssEnv$result_folder <-  result_folder
  ssEnv$result_folderData <-  dir_check_and_create(result_folder, "Data")
  ssEnv$result_folderChart <-    dir_check_and_create(result_folder, "Chart")
  ssEnv$result_folderInference <-    dir_check_and_create(result_folder, "Inference")
  ssEnv$result_folderPathway <-    dir_check_and_create(result_folder, "Pathway")
  ssEnv$result_folderPhenotype <-    dir_check_and_create(result_folder, "Phenotype")
  ssEnv$result_folderEuristic <-  dir_check_and_create(result_folder,"Euristic")
  ssEnv$session_folder <-  dir_check_and_create(result_folder,c("Log"))
  random_file_name <- paste(stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),".log", sep="")

  # Skip sink when running inside a callr child process — callr already
  # redirects stdout/stderr to log.txt; opening a second sink fights with
  # the redirect and can cause silent crashes on large I/O.
  if (!identical(Sys.getenv("SEMSEEKER_CHILD"), "1")) {
    if (sink.number() != 0)
      sink(NULL)
    file_name <- paste(as.character(Sys.info()["nodename"]),"_session_output.log", sep="")
    sink(file.path(ssEnv$session_folder,file_name), split = TRUE, append = TRUE)
  }

  foreachIndex <- 0

  # check if the arguments are valid
  arguments <- keys_create(ssEnv, arguments)
  ssEnv <- get_session_info()

  ssEnv$functionToExport <- c( "analyze_single_sample","deltar_single_sample",
    "dump_sample_as_bed_file", "delta_single_sample","dir_check_and_create",
    "file_path_build","analyze_single_sample_both",
    "sort_by_chr_and_start", "test_match_order", "lesions_get",
    "mutations_get")


  # to manage progress bar
  if(ssEnv$showprogress)
  {
    handler_settings <- progressr::handlers()
    if (!(exists("cli", mode = "function", inherits = TRUE)))
    {
      # check if handler is already registered
      if (!testthat::is_testing())
        # if (length(handler_settings$handler) == 0)
        if(!("cli" %in% handler_settings$handler))
        {
          progressr::handlers(global = TRUE)
          progressr::handlers("cli")
        }
    }
  }

  arguments <- set_env_variable(arguments,"maxResources",maxResources)
  arguments <- set_env_variable(arguments,"parallel_strategy","sequential")
  parallel_session()
  ssEnv <- get_session_info()

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will focus on:", paste(unique(ssEnv$keys_markers_figures$MARKER), collapse = " ", sep =" "),
    " due to ",  paste(unique(ssEnv$keys_markers_figures$FIGURE), collapse = " ", sep =" "),
    " of ",  paste( unique(ssEnv$keys_areas_subareas_markers_figures$AREA) , collapse = " ", sep =" "))

  # remove empty arguments
  v <- c()
  if(length(arguments)!=0)
  {
    temp <- arguments
    for (i in seq_along(temp))
    {
      if(is.null(temp[[i]]))
        v <- c(v,i)
      if(identical(temp[[i]],character(0)))
        v <- c(v,i)
    }
    # remove items with index in v
    if(length(v)!=0)
      arguments <- temp[-v]
    else
      arguments <- temp
  }
  # check length of arguments
  if(length(arguments)!=0)
  {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " This options are not recognized: ", paste(arguments, collapse = " ", sep =" "))
    # throw error
    stop("ERROR: This options are not recognized: ", paste(arguments, collapse = " ", sep =" "))
  }


  if(dry_run)
  {
    # out at console as pretty table

    knitr::kable(as.data.frame(ssEnv$keys_areas_subareas_markers_figures), format = "pipe", caption = "Selection:")
    message(ssEnv$keys_areas_subareas_markers_figures)
    stop("INFO: Dry run is requested. Exiting now.")
  }

  update_session_info(ssEnv)

  return(ssEnv)
}
