# Standalone debug script — NOT a testthat file.
# Goal: verify depth=3 with kruskal.test actually produces DEPTH>1 rows,
# bypassing the slow full test file.
# Run from package root: Rscript tests/aa_debug_t1_kruskal.R

library(semseeker)
source("tests/testthat/setup.R")

# ── Inline copy of .aa_setup_result_folder ────────────────────────────────────
seed <- 777
n_probes <- 200L
areas <- c("GENE", "POSITION")

tempFolder <- tempFolders[1]
tempFolders <- tempFolders[-1]
unlink(tempFolder, recursive = TRUE)

set.seed(seed)
local_probes <- probe_features[seq_len(n_probes), ]
local_sig <- matrix(stats::rbeta(n_probes * nsamples, 90L, 10L),
                     nrow = n_probes, ncol = nsamples)
local_sig[1:50, 1:5] <- stats::rbeta(50L * 5L, 1L, 100L)
rownames(local_sig) <- local_probes$PROBE
local_sig <- as.data.frame(local_sig)
colnames(local_sig) <- mySampleSheet$Sample_ID

cat("\n=== STEP 1: semseeker (areas=GENE+POSITION, markers=MUTATIONS) ===\n")
t0 <- Sys.time()
semseeker::semseeker(
  input             = local_sig,
  sample_sheet      = mySampleSheet,
  result_folder     = tempFolder,
  parallel_strategy = parallel_strategy,
  areas             = areas,
  markers           = c("MUTATIONS"),
  start_fresh       = TRUE,
  showprogress      = FALSE,
  verbosity         = 1
)
cat(sprintf("semseeker took %.1f sec\n", as.numeric(Sys.time() - t0, units = "secs")))

# ── INSPECT: what pivot parquet files were created? ──────────────────────────
cat("\n=== STEP 1b: list parquet pivot files produced ===\n")
parquets <- list.files(tempFolder, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
cat(length(parquets), "parquet files:\n")
print(basename(parquets))

# ── STEP 2: association_analysis with KW + depth=3 ───────────────────────────
inference_details <- data.frame(
  independent_variable = "Sample_Group",      # categorical (3 groups)
  family_test          = "kruskal.test",      # exact string from R/test_model.R:93
  transformation_y     = "",
  transformation_x     = "",
  depth_analysis       = 3L,
  filter_p_value       = FALSE,
  stringsAsFactors     = FALSE
)

cat("\n=== STEP 2: association_analysis(depth=3, kruskal.test) ===\n")
t1 <- Sys.time()
semseeker:::association_analysis(
  inference_details = inference_details,
  result_folder     = tempFolder,
  parallel_strategy = parallel_strategy,
  markers           = c("MUTATIONS"),
  figures           = c("HYPO"),
  multiple_test_adj = "BH",
  showprogress      = FALSE,
  verbosity         = 1
)
cat(sprintf("association_analysis took %.1f sec\n", as.numeric(Sys.time() - t1, units = "secs")))

# ── STEP 3: inspect output CSVs ──────────────────────────────────────────────
cat("\n=== STEP 3: result CSVs and DEPTH distribution ===\n")
inference_dir <- file.path(tempFolder, "Inference")
csv_files <- list.files(inference_dir, pattern = "\\.csv$", recursive = TRUE,
                         full.names = TRUE)
cat(length(csv_files), "CSV files:\n")
for (f in csv_files) {
  cat(sprintf("  - %s  (%d bytes)\n", basename(f), file.info(f)$size))
}

result_csv <- csv_files[!grepl("(?i)covariates_model", csv_files)][1]
if (!is.na(result_csv) && file.exists(result_csv) && file.info(result_csv)$size > 10) {
  df <- utils::read.csv2(result_csv)
  cat(sprintf("\nMain result CSV: %s — %d rows, %d cols\n", basename(result_csv), nrow(df), ncol(df)))
  if ("DEPTH" %in% colnames(df)) {
    cat("DEPTH value distribution:\n")
    print(table(df$DEPTH))
  } else {
    cat("(no DEPTH column in output)\n")
  }
}

cat("\n=== DONE ===\n")
unlink(tempFolder, recursive = TRUE)
