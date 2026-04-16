parallel_session <- function()
{
  #
  ssEnv <- get_session_info()
  parallel_strategy <- ssEnv$parallel_strategy

  # ---- Guard: fork() + Polars C++ thread pool on macOS → silent crash -----
  # `multicore` uses fork() on macOS. If Polars has already been loaded in
  # the parent process (e.g. by a prior signal_save() call), its C++ thread
  # pool makes fork() undefined behaviour — the child is killed by a Mach
  # exception with no R-visible error. Force multisession to avoid this.
  if (parallel_strategy == "multicore" &&
      identical(Sys.info()[["sysname"]], "Darwin") &&
      "polars" %in% loadedNamespaces())
  {
    log_event("WARNING: parallel_strategy='multicore' + Polars on macOS is ",
              "unsafe (fork() on a process with active C++ thread pool ",
              "causes silent crashes). Forcing parallel_strategy='multisession'.")
    parallel_strategy <- "multisession"
    ssEnv$parallel_strategy <- "multisession"
  }

  # check if os is macos
  if (Sys.info()["sysname"] == "Darwin" & parallel_strategy != "sequential")
  {
    env_var <- Sys.getenv("OBJC_DISABLE_INITIALIZE_FORK_SAFETY")
    if (env_var != "YES")
    {
      log_event("ERROR: Setting OBJC_DISABLE_INITIALIZE_FORK_SAFETY must be YES to work in multiprocess. \n
        execute: export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES at shell or
         Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES') ! ")
      # Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY="YES")
      stop("Setting OBJC_DISABLE_INITIALIZE_FORK_SAFETY must be YES to work in multiprocess. \n
        execute: export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES at shell or
         Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES') !")
    }
  }

  if(parallelly::supportsMulticore())
    options(parallelly.fork.enable= TRUE)
  else
    options(parallelly.fork.enable= FALSE)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    nCore <- 2L
  } else {
    # use all cores in devtools::test()
    nCore <- future::availableCores() - 1
    nCore <- if(floor(future::availableCores() * ssEnv$maxResources/100 ) > nCore ) nCore else floor(future::availableCores() * ssEnv$maxResources/100 )
    nCore <- max(1L, nCore)  # guarantee at least 1 worker (e.g. covr subprocess with 1 core)
  }
  # permutation cluster
  outFile <- file.path(ssEnv$session_folder, "cluster_r.out")

  #
  ssEnv$parallel <- data.frame("parallel_strategy"="")
  ssEnv$parallel$parallel_strategy <- parallel_strategy
  ssEnv$parallel$nCore <- nCore

  options(readr.num_threads = nCore)  # Set to the desired number of threads

  options(doFuture.foreach.export = ".export-and-automatic-with-warning")
  # doFuture
  # is changed in onload see zzz.R
  # options(parallelly.fork.enable= FALSE)
  # options(future.globals.resolve = TRUE)
  # allow export of object of 32gb with future
  options(future.globals.maxSize= 64 * 1024^3)


  # check if future is registered
  # backend_name <- foreach::getDoParName()
  # if (!is.null(backend_name)) {
  #   log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " doFuture is registered as the %dopar% backend.")
  # } else {
  #   log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " doFuture is not registered as the %dopar% backend.")
  #   doFuture::registerDoFuture()
  # }
  doFuture::registerDoFuture()

  # get the future plan
  future_plan <- future::plan()
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Future Plan: ", future_plan)

  # TODO: improve planning parallel management using also cluster
  if(parallel_strategy=="multisession")
  {
    future::plan( future::multisession, workers = nCore)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in multisession with:", nCore, " Cores")
  }
  if(parallel_strategy=="multicore")
  {
    future::plan( future::multicore, workers = nCore)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in multicore with:", nCore," Cores")
  }
  if(parallel_strategy=="cluster")
  {
    if (!is.null(ssEnv$cluster_workers))
    {
      future::plan( future::cluster, workers = ssEnv$cluster_workers)
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:",ssEnv$cluster_workers)
    }
    else
    {
      future::plan( future::cluster, workers = nCore)
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work with a cluster with:", nCore," Cores")
    }
  }

  if(parallel_strategy!="multisession" & parallel_strategy!="multicore"
    & parallel_strategy!="cluster")
  {
    options(parallelly.fork.enable= FALSE)
    future::plan(strategy = future::sequential)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work in sequential mode")
  }

  update_session_info(ssEnv)

}
