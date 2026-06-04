parallel_session <- function()
{
  #
  ssEnv <- get_session_info()
  parallel_strategy <- ssEnv$parallel_strategy

  # NOTE: `multicore` on macOS uses fork() and is known to be unsafe in
  # combination with Polars' C++ thread pool — forked children can be
  # killed by a Mach exception with no R-visible error. Tests on macOS
  # now default to `multisession` (see setup.R). End users on macOS
  # should use `multisession` or `sequential`.
  #
  # E-14: `multisession` workers are fresh R processes — .pkgglobalenv$ssEnv
  # starts empty. Every %dorng% foreach body must call
  # update_session_info(ssEnv) as its first statement to populate the
  # worker's namespace. See engineering-decisions.md §1.3.

  # macOS: OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES is recommended for any
  # non-sequential strategy (some packages still attempt fork internally).
  # Bioconductor disallows Sys.setenv() in package code; users on macOS
  # using a non-sequential strategy must set the env var themselves before
  # calling init_env(), e.g. in ~/.Renviron or Sys.setenv() at the top of
  # their script. We log a warning when the env var is missing so the
  # cause of any subsequent fork crash is obvious.
  if (Sys.info()["sysname"] == "Darwin" && parallel_strategy != "sequential") {
    env_var <- Sys.getenv("OBJC_DISABLE_INITIALIZE_FORK_SAFETY")
    if (env_var != "YES") {
      log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " OBJC_DISABLE_INITIALIZE_FORK_SAFETY is not set to YES.",
                " On macOS with a non-sequential parallel strategy this can",
                " cause fork-related crashes. Set it in ~/.Renviron or call",
                " Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = 'YES')",
                " before init_env().")
    }
  }

  # E-13: multicore (fork) on macOS is unsafe with Polars' C++ thread pool —
  # forked children are killed by Mach exceptions with no R-visible error.
  # Force multisession (separate R processes) instead.
  if (Sys.info()["sysname"] == "Darwin" && parallel_strategy == "multicore") {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " multicore (fork) is unsafe on macOS with Polars. Switching to multisession.")
    parallel_strategy <- "multisession"
    ssEnv$parallel_strategy <- parallel_strategy
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
