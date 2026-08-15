# ---------------------------------------------------------------------------
# Windows VM — one-off setup, then the fast loop
#
# Why this exists. The Windows CI job takes 4-5 hours and, when it dies, can
# leave nothing behind: run 31852499666 failed with no failed step and no log
# (`BlobNotFound`). A local Windows VM does not make a full R CMD check faster —
# on Apple Silicon R runs under x64 emulation, so it is slower — but it gives
# the two things the CI cannot: the output, and the ability to re-run ONE test
# file in minutes instead of waiting for a whole matrix.
#
# Host: Apple Silicon -> Parallels or VMware Fusion with Windows 11 ARM.
# VirtualBox does not run Windows on ARM. R for Windows is x64 and runs under
# the OS emulation layer; for the defects Windows actually produces — path
# separators, case-insensitive file systems, file locking, encoding, `parallel`
# without fork — that is equivalent to the x64 runner. For a genuinely
# architecture-specific defect it is not, and only the CI can tell you.
#
# Inside the VM, with the repo reachable (Parallels shares the Mac home under
# \\Mac\Home), run this file once with:
#
#   "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" dev\windows-vm-setup.R
#
# ---------------------------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

message("R: ", R.version.string, " on ", R.version$platform)
if (.Platform$OS.type != "windows")
  stop("This script is meant to be run inside the Windows VM.")

## 1. Toolchain -------------------------------------------------------------
## Rtools is required to build the few packages that have no Windows binary.
if (!nzchar(Sys.which("make")))
  message("NOTE: Rtools does not seem to be on PATH. Install it from ",
          "https://cran.r-project.org/bin/windows/Rtools/ and reopen the shell.")

## 2. Package managers ------------------------------------------------------
for (p in c("remotes", "BiocManager", "devtools", "testthat"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

## 3. Dependencies ----------------------------------------------------------
## The long pole, but it is paid ONCE: on Windows these are mostly binaries.
## BiocManager is used for everything so the Bioconductor annotation packages
## resolve; limma goes first for the reason recorded in the CI workflow — the
## Illumina annotation packages load minfi while installing, and minfi needs it.
BiocManager::install("limma", ask = FALSE, update = FALSE)

deps <- remotes::local_package_deps(dependencies = TRUE)
missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing ", length(missing), " dependencies: ",
          paste(utils::head(missing, 10), collapse = ", "),
          if (length(missing) > 10) ", ..." else "")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

## polars is not on CRAN and needs NOT_CRAN, same as the workflow does.
Sys.setenv(NOT_CRAN = "true")
if (!requireNamespace("polars", quietly = TRUE))
  install.packages("polars",
                   repos = "https://community.r-multiverse.org")

message("\nSetup done. The loop from here on is:\n\n",
        "  devtools::install('.', quick = TRUE, upgrade = FALSE, build = FALSE)\n",
        "  testthat::test_dir('tests/testthat', package = 'SEMseeker',\n",
        "                     load_package = 'installed', stop_on_failure = FALSE)\n\n",
        "and to iterate on one file, add filter = '7-association_analysis'.\n\n",
        "Use load_package = 'installed' and not devtools::test(): under load_all\n",
        "a helper's closure is package:SEMseeker, whose environment chain reaches\n",
        "neither the global environment nor setup.R's, so the two harnesses do not\n",
        "measure the same thing (see AI-254).\n")
