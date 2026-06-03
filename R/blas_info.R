# AI-060: BLAS detection + setup helper.
#
# Two things this file does:
#   (B) .blas_info() introspects extSoftVersion()['BLAS'] and decides
#       whether the BLAS R is linked against is multi-threaded. init_env
#       emits a one-line WARNING when it isn't — the AI-040 limma/voom
#       path spends most of its wall-clock in solve()/crossprod() and
#       loses ~5-10x on machines using the single-thread reference BLAS.
#   (C) setup_blas() is a user-facing helper that prints the current
#       BLAS config and, on macOS, the exact one-time command to switch
#       to Accelerate (vecLib). Read-only by default — never modifies
#       the system. If apply = TRUE it also tries to apply the fix via
#       a sudo prompt.
#
# Design: no new package dependency. Detection is purely from
# extSoftVersion()['BLAS'] string match — fast, no RhpcBLASctl needed.

.blas_info <- function() {
  blas_path <- tryCatch(unname(extSoftVersion()["BLAS"]),
                         error = function(e) "")
  if (!is.character(blas_path) || length(blas_path) == 0L) blas_path <- ""

  is_accelerate <- grepl("vecLib|Accelerate", blas_path, ignore.case = TRUE)
  is_openblas   <- grepl("openblas",          blas_path, ignore.case = TRUE)
  is_mkl        <- grepl("mkl",               blas_path, ignore.case = TRUE)
  is_reference  <- grepl("libRblas(\\.[0-9])?\\.(dylib|so)$", blas_path)

  flavor <- if (is_accelerate) "Accelerate (vecLib)"
            else if (is_openblas)  "OpenBLAS"
            else if (is_mkl)       "MKL"
            else if (is_reference) "reference (single-thread)"
            else "unknown"

  list(path           = blas_path,
       multi_threaded = is_accelerate || is_openblas || is_mkl,
       flavor         = flavor,
       is_reference   = is_reference)
}

.warn_blas_single_thread <- function() {
  info <- .blas_info()
  if (!info$multi_threaded) {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " BLAS is ", info$flavor,
              " (", info$path, "). limma/voom families will run",
              " single-thread — typically 5-10x slower than with",
              " a multi-threaded BLAS. Run SEMseeker::setup_blas()",
              " for instructions.")
  }
  invisible(info)
}

#' Print BLAS configuration and (optionally) apply the macOS speedup
#'
#' SEMseeker's batch families (\code{limma_<N>}, \code{voom_<N>}, added
#' in AI-040) spend most of their wall-clock inside BLAS calls
#' (\code{solve()}, \code{crossprod()} on the design matrix). With a
#' multi-threaded BLAS those calls use every physical core; with the
#' single-thread reference BLAS R ships by default on macOS they run
#' on one core. The difference is typically 5-10x on Apple Silicon.
#'
#' Read-only by default. Call \code{setup_blas(apply = TRUE)} to print
#' the macOS fix command and ATTEMPT to run it with sudo (requires
#' your password).
#'
#' @param apply logical. If TRUE and the platform is macOS with the
#'   reference BLAS, prompts for sudo and runs the symlink command.
#' @return Invisible list with current BLAS info: \code{path},
#'   \code{flavor}, \code{multi_threaded}, \code{is_reference}.
#'
#' @examples
#' \dontrun{
#'   SEMseeker::setup_blas()                 # just print info
#'   SEMseeker::setup_blas(apply = TRUE)     # try to fix on macOS
#' }
#' @export
setup_blas <- function(apply = FALSE) {
  info <- .blas_info()
  cat("BLAS path:      ", info$path,            "\n", sep = "")
  cat("Flavor:         ", info$flavor,          "\n", sep = "")
  cat("Multi-threaded: ", info$multi_threaded,  "\n", sep = "")
  if (info$multi_threaded) {
    cat("\nNothing to do — BLAS is already multi-threaded.\n")
    return(invisible(info))
  }

  if (Sys.info()[["sysname"]] != "Darwin") {
    cat("\nNon-macOS system. The vecLib trick below applies only to macOS.\n",
        "Linux: install OpenBLAS or MKL via your package manager and\n",
        "       relink R (Rcmd config BLAS_LIBS).\n",
        "Windows: install Microsoft R Open or use Rtools with MKL.\n",
        sep = "")
    return(invisible(info))
  }

  rhome <- R.home()
  cmd <- sprintf(
    "cd %s/lib && sudo ln -sf libRblas.vecLib.dylib libRblas.dylib",
    rhome)
  cat("\nmacOS: ship-with R uses the single-thread reference BLAS by default.\n",
      "Switching to Accelerate (vecLib) typically gives 5-10x speedup for\n",
      "limma::lmFit and voom on >100k areas. One-time fix:\n\n  ",
      cmd, "\n\nThen restart R. Verify with extSoftVersion()['BLAS'] — it\n",
      "should contain 'vecLib'.\n",
      sep = "")

  if (isTRUE(apply)) {
    cat("\nAttempting to run the symlink now (will prompt for sudo) …\n")
    rc <- tryCatch(system(cmd, intern = FALSE), error = function(e) -1L)
    if (rc == 0L) {
      cat("Symlink applied. Restart R, then check setup_blas() again.\n")
    } else {
      cat("Symlink command exited with status ", rc,
          " — please run it manually.\n", sep = "")
    }
  }
  invisible(info)
}
