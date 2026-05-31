# Diagnostic — why does depth=3 silently write only DEPTH=1 rows?
# Hypothesis to test: pivot_file_name_parquet() builds a path that doesn't match
# the actual file on disk OR ssEnv$keys_areas_subareas_markers_figures is empty
# for MARKER=MUTATIONS.
# Run: Rscript tests/aa_debug_depth3.R

library(semseeker)
source("tests/testthat/setup.R")

# Reduced to 50 probes for speed
n_probes <- 50L
tempFolder <- tempFolders[1]; tempFolders <- tempFolders[-1]
unlink(tempFolder, recursive = TRUE)

set.seed(777)
local_probes <- probe_features[seq_len(n_probes), ]
local_sig <- matrix(stats::rbeta(n_probes * nsamples, 90L, 10L),
                     nrow = n_probes, ncol = nsamples)
local_sig[1:25, 1:5] <- stats::rbeta(25L * 5L, 1L, 100L)
rownames(local_sig) <- local_probes$PROBE
local_sig <- as.data.frame(local_sig)
colnames(local_sig) <- mySampleSheet$Sample_ID

cat("\n=== STEP 1: semseeker ===\n")
t0 <- Sys.time()
semseeker::semseeker(
  input             = local_sig,
  sample_sheet      = mySampleSheet,
  result_folder     = tempFolder,
  parallel_strategy = parallel_strategy,
  areas             = c("GENE", "POSITION"),
  markers           = c("MUTATIONS"),
  start_fresh       = TRUE,
  showprogress      = FALSE,
  verbosity         = 1
)
cat(sprintf("semseeker: %.1f sec\n", as.numeric(Sys.time() - t0, units = "secs")))

# ── STEP 2: inspect ssEnv state AFTER semseeker (must reinit) ────────────────
cat("\n=== STEP 2: re-init ssEnv (semseeker closes it on exit) ===\n")
semseeker:::init_env(tempFolder, parallel_strategy = "sequential",
                     showprogress = FALSE, verbosity = 1)
ssEnv <- semseeker:::get_session_info()
cat("ssEnv$genome_build =", paste(ssEnv$genome_build), "\n")
cat("ssEnv$tech         =", paste(ssEnv$tech), "\n")
cat("ssEnv has keys_areas_subareas_markers_figures? ",
    "keys_areas_subareas_markers_figures" %in% names(ssEnv), "\n")
keys_all <- ssEnv$keys_areas_subareas_markers_figures
if (!is.null(keys_all)) {
  cat("Total keys rows:", nrow(keys_all), "\n")
  cat("MARKER counts:\n"); print(table(keys_all$MARKER))
  mut <- keys_all[keys_all$MARKER == "MUTATIONS", ]
  cat("\nMUTATIONS subset rows:", nrow(mut), "\n")
  cat("MUTATIONS AREA x SUBAREA combos:\n")
  print(unique(mut[, c("FIGURE","AREA","SUBAREA")]))
} else {
  cat("!! keys_areas_subareas_markers_figures is NULL — depth>1 cannot run\n")
}

# ── STEP 3: actual disk vs constructed-path match ────────────────────────────
cat("\n=== STEP 3: file.exists() per MUTATIONS GENE key ===\n")
if (!is.null(keys_all)) {
  mut_gene <- keys_all[keys_all$MARKER == "MUTATIONS" & keys_all$AREA == "GENE", ]
  for (k in seq_len(min(nrow(mut_gene), 5))) {
    row <- mut_gene[k, ]
    path <- semseeker:::pivot_file_name_parquet(row$MARKER, row$FIGURE, row$AREA, row$SUBAREA)
    cat(sprintf("  %s_%s_%s_%s\n    expected:%s\n    exists  :%s\n",
                row$MARKER, row$FIGURE, row$AREA, row$SUBAREA,
                basename(path), file.exists(path)))
  }
}

# ── STEP 4: actual parquet files on disk ─────────────────────────────────────
cat("\n=== STEP 4: parquet files on disk ===\n")
parquets <- list.files(tempFolder, pattern = "\\.parquet$",
                        recursive = TRUE, full.names = FALSE)
cat(length(parquets), "files; head 6:\n")
print(head(parquets, 6))

# ── STEP 5: read the session log to see WARNINGs that fired ──────────────────
cat("\n=== STEP 5: tail of session log (looking for WARNING/depth) ===\n")
log_files <- list.files(tempFolder, pattern = "session_output\\.log$",
                         recursive = TRUE, full.names = TRUE)
cat("log files:", length(log_files), "\n")
for (lf in log_files) {
  cat("\n--- ", lf, " (last 30 WARNING/depth lines) ---\n")
  lines <- readLines(lf, warn = FALSE)
  rel <- grep("WARNING|depth|pivot|File not found|skipped", lines, ignore.case = TRUE, value = TRUE)
  cat(paste(tail(rel, 30), collapse = "\n"), "\n")
}

cat("\n=== DONE — tempFolder retained for inspection:\n", tempFolder, "\n")
