# num_rows <- 3e^6
# num_cols <- 5200
# populationMatrix <- as.data.frame(matrix(runif(num_rows * num_cols), nrow = num_rows, ncol = num_cols))
# library(future)
# options(future.globals.maxSize = 10 * 1024^3)
Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES')

rm(list = ls())
# DEBUG: trace setup.R progress so we can see (in the macOS run log) at
# which step GEOquery / tcltk loading is triggered. AI-017 tracking.
.trace_step <- function(msg) {
  cat(sprintf("[SETUP-TRACE %s] %s\n",
              format(Sys.time(), "%H:%M:%OS3"), msg))
  flush.console()
}
.trace_step("setup.R BEGIN")
loadNamespace("future")
.trace_step("after loadNamespace(future)")
loadNamespace("stats")
.trace_step("after loadNamespace(stats)")
Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = 'YES')

# ── Probe features: 20k real EPIC probe IDs (imprinting DMR-aware) ────────
.trace_step("loading bundled test_master_features fixture")
utils::data("test_master_features", package = "SEMseeker", envir = environment())
.probe_features_all <- as.data.frame(test_master_features, stringsAsFactors = FALSE)
.trace_step(sprintf("loaded test_master_features: %d probes", nrow(.probe_features_all)))

# ── Beta matrix: real GSE133774 (EPIC 850k, BWS + MLID, 10 samples) ───────
# Built once by data-raw/build_test_signal_fixture.R; ships in data/.
# Single source of truth for both automated tests and vignette (AI-123).
.trace_step("loading bundled test_signal_gse133774 fixture")
utils::data("test_signal_gse133774",      package = "SEMseeker", envir = environment())
utils::data("test_samplesheet_gse133774", package = "SEMseeker", envir = environment())

# Align probe_features to probes present in the beta matrix
.common <- intersect(.probe_features_all$PROBE, rownames(test_signal_gse133774))
probe_features <- .probe_features_all[.probe_features_all$PROBE %in% .common, ]
signal_data    <- as.data.frame(test_signal_gse133774[.common, , drop = FALSE])

nprobes  <<- nrow(signal_data)
nsamples <<- ncol(signal_data)
.trace_step(sprintf("aligned fixture: %d probes × %d samples", nprobes, nsamples))

# ── Sample sheet ───────────────────────────────────────────────────────────
mySampleSheet <- test_samplesheet_gse133774
# Extra columns used by association / batch tests
set.seed(474693)
mySampleSheet$Phenotest   <- stats::rnorm(nrow(mySampleSheet), mean = 1000, sd = 567)
mySampleSheet$Group       <- c(rep(TRUE,  ceiling(nrow(mySampleSheet) / 2)),
                               rep(FALSE, floor(nrow(mySampleSheet)   / 2)))
mySampleSheet$Covariates1 <- stats::rnorm(nrow(mySampleSheet), mean = 567,  sd = 1000)
mySampleSheet$Covariates2 <- stats::rnorm(nrow(mySampleSheet), mean = 67,   sd = 100)

mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
signal_data_batch   <<- list(signal_data,   signal_data,   signal_data)

# ── IQR thresholds from Reference samples only ────────────────────────────
.ref_ids     <- unique(mySampleSheet$Sample_ID[mySampleSheet$Sample_Group == "Reference"])
.ref_ids     <- .ref_ids[.ref_ids %in% colnames(signal_data)]
.ref_signal  <- signal_data[, .ref_ids, drop = FALSE]
q1           <- apply(.ref_signal, 1, function(x) stats::quantile(x, 0.25, na.rm = TRUE))
q3           <- apply(.ref_signal, 1, function(x) stats::quantile(x, 0.75, na.rm = TRUE))
signal_medians <- apply(.ref_signal, 1, stats::median)
iqr          <- data.frame(q3 - q1)

signal_superior_thresholds <- data.frame("HIGH" = q3 + 3 * iqr)
signal_inferior_thresholds <- data.frame("LOW"  = q1 - 3 * iqr)
colnames(signal_inferior_thresholds) <- "LOW"
colnames(signal_superior_thresholds) <- "HIGH"
row.names(signal_superior_thresholds) <- probe_features$PROBE
row.names(signal_inferior_thresholds) <- probe_features$PROBE

signal_thresholds <- data.frame(
  "signal_median_values"       = signal_medians,
  "signal_inferior_thresholds" = signal_inferior_thresholds,
  "signal_superior_thresholds" = signal_superior_thresholds,
  "iqr"   = iqr,
  "q1"    = q1,
  "q3"    = q3
)
colnames(signal_thresholds) <- c("signal_median_values", "signal_inferior_thresholds",
                                  "signal_superior_thresholds", "iqr", "q1", "q3")
signal_thresholds$CHR   <- probe_features$CHR
signal_thresholds$START <- probe_features$START
signal_thresholds$END   <- probe_features$END

mySampleSheet              <<- mySampleSheet
signal_data                <<- signal_data
signal_medians             <<- signal_medians
signal_inferior_thresholds <<- signal_inferior_thresholds
signal_superior_thresholds <<- signal_superior_thresholds
signal_thresholds          <<- signal_thresholds
nsamples                   <<- nsamples
nprobes                    <<- nprobes



LESIONS_BP <<- 5000L  # AI-092 + AI-044 merged: bp-based window, literature-aligned default 5 kbp (see AI-048).
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3
# "multicore" (fork) is unsafe on macOS with Polars' C++ thread pool — forked
# children can be killed by Mach exceptions.  Use "multisession" when the
# package is installed (R CMD check, CI, devtools::install()).  Fall back to
# "multicore" only when running under devtools::load_all() on non-macOS,
# because multisession workers cannot see load_all()'d internals.
if (Sys.info()[["sysname"]] == "Darwin") {
  parallel_strategy <<- "multisession"
} else {
  parallel_strategy <<- "multicore"
}
markers <<- c("MUTATIONS","DELTAQ","DELTARQ","DELTAP","DELTARP","LESIONS")


# TODO
# recover session stored
#
tmp <- normalizePath(tempdir())
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")


check_execution_context <- function() {
  calls <- sys.calls()
  if (any(sapply(calls, function(x) "test_file" %in% names(x)))) {
    showprogress <<- FALSE
    verbosity <<- 1
    # message("Called from testthat")
  } else {
    showprogress <<- TRUE
    verbosity <<- 4
    # message("Called from source or directly")
  }
}

check_execution_context()

# TODO
# recover session stored
#
tmp <- tempdir()
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")
