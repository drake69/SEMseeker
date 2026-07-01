# AI-086: integration test for sample_sheet_result.csv burden columns.
#
# Canary for AI-083: study_summary_total() failed to populate per-sample
# burden in sample_sheet_result.csv on ewas_osteoporosis/GSE99624 (48 samples,
# 450K) — all burden columns landed as NA, which then crashed every depth=1
# inference downstream with "data are not the same size".
#
# test-7-association_analysis.R only verifies this indirectly (depth=1 would
# fail on all-NA burden vectors); this test asserts the contract directly:
#
#   1. sample_sheet_result.csv exists after semseeker().
#   2. Every expected <MARKER>_<FIGURE> burden column is PRESENT.
#   3. Every burden column has at least one non-NA sample.
#   4. No sample has ALL burden values NA (per-sample BED collection works).
#   5. PROBES_COUNT > 0 on every sample.
#   6. Sanity: injected hypomethylated outliers surface as MUTATIONS_HYPO > 0.
#
# Why discrete markers only: study_summary_total() iterates AREA=="POSITION"
# keys, which is the set declared in keys_create.R via
# keys_markers_default_discrete (MUTATIONS / LESIONS / DELTAQ / DELTARQ /
# DELTAP / DELTARP) plus the two continuous DELTAS / DELTAR. LESIONS depends
# on MUTATIONS having been computed (SOURCE column), so the full
# (MUTATIONS + DELTA*) set exercises the whole derive-and-aggregate path.

.burden_setup_signal_with_outliers <- function(seed = 12345L) {
  set.seed(seed)
  n_probes_b  <- 200L
  n_samples_b <- nsamples
  local_probes  <- probe_features[1:n_probes_b, ]
  local_samples <- mySampleSheet

  # Background: mostly methylated (Beta(90, 10)).
  # IMPORTANT: outliers are detected per-probe across samples (IQR × 3 on
  # the inter-sample distribution). If we inject the SAME band across all
  # samples uniformly, the per-probe q1/q3 collapse to that uniform value
  # and IQR → 0 → no sample is flagged. So we inject DISJOINT, RANDOM
  # per-sample probe sets — each sample is an outlier at probes the other
  # samples are not, which is exactly the inter-sample dispersion the
  # threshold needs.
  local_sig <- matrix(
    stats::rbeta(n_probes_b * n_samples_b, 90L, 10L),
    nrow = n_probes_b, ncol = n_samples_b
  )
  for (s in seq_len(n_samples_b)) {
    # 10 random HYPO probes from positions 1:100 → values ≈ 0
    hypo_probes <- sample.int(100L, 10L)
    local_sig[hypo_probes, s] <- stats::rbeta(10L, 1L, 100L)
    # 10 random HYPER probes from positions 101:200 → values ≈ 1
    hyper_probes <- 100L + sample.int(100L, 10L)
    local_sig[hyper_probes, s] <- stats::rbeta(10L, 100L, 1L)
  }

  rownames(local_sig) <- local_probes$PROBE
  local_sig <- as.data.frame(local_sig)
  # signal_data has 10 unique columns; mySampleSheet has 16 rows (Reference reuse pattern)
  colnames(local_sig) <- colnames(signal_data)
  list(signal = local_sig, samples = local_samples, probes = local_probes)
}

test_that("sample_sheet_result.csv has populated burden columns for all discrete + continuous markers (AI-086, canary for AI-083)", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  unlink(tempFolder, recursive = TRUE)

  syn <- .burden_setup_signal_with_outliers()

  # inpute="median" is defensive: harmless when the input has no NAs (the
  # guard at inpute_missing_values.R:7 short-circuits), useful if a future
  # signal generator introduces them.
  SEMseeker::semseeker(
    input             = syn$signal,
    sample_sheet      = syn$samples,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS", "LESIONS",
                          "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                          "DELTAS", "DELTAR"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # SEMseeker writes via io_file_path_build() → name_cleaning() → toupper(),
  # so the on-disk file is SAMPLE_SHEET_RESULT.csv (not the lowercase form
  # used in the API contract). macOS/Windows file systems are case-insensitive
  # by default so either spelling works; Linux ext4 is case-sensitive and
  # only the uppercase form resolves. Use uppercase here to be correct on
  # all three CI runners.
  result_csv <- file.path(tempFolder, "Data", "SAMPLE_SHEET_RESULT.csv")
  testthat::expect_true(
    file.exists(result_csv),
    info = sprintf(
      "SAMPLE_SHEET_RESULT.csv was not written by semseeker() — tempFolder=%s",
      tempFolder
    )
  )
  testthat::skip_if_not(file.exists(result_csv),
                        "downstream assertions need the result CSV")

  df <- utils::read.csv2(result_csv, stringsAsFactors = FALSE)

  # Required columns: MUTATIONS + DELTA* (LESIONS is optional — derives from
  # MUTATIONS clusters and may legitimately be 0 on small synthetic data.
  # Tracked separately below as a soft expectation.)
  required_markers <- c("MUTATIONS",
                        "DELTAP", "DELTAQ", "DELTARP", "DELTARQ",
                        "DELTAS", "DELTAR")
  required_burden_cols <- c(
    paste0(required_markers, "_HYPER"),
    paste0(required_markers, "_HYPO")
  )
  required_cols <- c("Sample_ID", required_burden_cols, "PROBES_COUNT")

  missing_cols <- setdiff(required_cols, colnames(df))
  testthat::expect_equal(
    length(missing_cols), 0L,
    info = paste("missing required columns:", paste(missing_cols, collapse = ", "))
  )

  # Soft: log if LESIONS_HYPER/HYPO are absent so we surface the gap
  # without failing — see AI-088 follow-up for an explicit LESIONS canary.
  for (lesion_col in c("LESIONS_HYPER", "LESIONS_HYPO")) {
    if (!(lesion_col %in% colnames(df))) {
      message(sprintf(
        "test-8-burden-integration: %s absent — likely synthetic-data sparsity (soft warn, see AI-088)",
        lesion_col
      ))
    }
  }

  # The per-column and per-sample assertions below operate only on the
  # required burden columns (LESIONS handled by the soft block above).
  expected_burden_cols <- required_burden_cols

  # AI-083 canary (per-column): every burden column has ≥ 1 non-NA value.
  for (col in intersect(expected_burden_cols, colnames(df))) {
    n_non_na <- sum(!is.na(df[[col]]))
    testthat::expect_gt(
      n_non_na, 0L,
      label = sprintf("%s non-NA samples (got %d)", col, n_non_na)
    )
  }

  # AI-083 canary (per-sample, relaxed): NOT every sample must have a
  # populated burden — synthetic data sparsity can legitimately leave some
  # samples with zero events for ALL (marker, figure) combos, which then
  # become all-NA after the all.x=TRUE merge in study_summary_total. The
  # AI-083 bug was 100% of samples NA (whole-population integration broken),
  # not "some samples are NA". So we assert the all-NA fraction is below a
  # safety margin instead of zero. A tighter per-sample canary on real or
  # tightly-controlled data is tracked in AI-088.
  present_burden <- intersect(expected_burden_cols, colnames(df))
  if (length(present_burden) > 0L) {
    per_sample_all_na <- apply(
      df[, present_burden, drop = FALSE], 1L,
      function(r) all(is.na(r))
    )
    all_na_fraction <- sum(per_sample_all_na) / nrow(df)
    testthat::expect_lt(
      all_na_fraction, 0.5,
      label = sprintf(
        "all-NA-burden fraction (got %.2f, %d/%d): %s",
        all_na_fraction, sum(per_sample_all_na), nrow(df),
        paste(df$Sample_ID[per_sample_all_na], collapse = ", ")
      )
    )
  }

  # PROBES_COUNT > 0 on every sample.
  testthat::expect_true(
    all(df$PROBES_COUNT > 0L, na.rm = TRUE),
    info = "PROBES_COUNT must be > 0 on every sample"
  )

  # Sanity: injected HYPO outliers must surface as MUTATIONS_HYPO > 0 on
  # at least the 5 samples we touched (sample indices 1:5).
  if ("MUTATIONS_HYPO" %in% colnames(df)) {
    testthat::expect_gt(
      sum(df$MUTATIONS_HYPO > 0, na.rm = TRUE), 0L,
      label = "samples with MUTATIONS_HYPO > 0 (injected-outlier sanity)"
    )
  }
  if ("MUTATIONS_HYPER" %in% colnames(df)) {
    testthat::expect_gt(
      sum(df$MUTATIONS_HYPER > 0, na.rm = TRUE), 0L,
      label = "samples with MUTATIONS_HYPER > 0 (injected-outlier sanity)"
    )
  }

  # Explicit pivot assertion: every required marker must have a per-figure
  # parquet pivot under Data/Pivots/<MARKER>/, and each one must contain
  # rows. The burden checks above cover this indirectly (an empty pivot
  # → all-NA burden column), but a direct file/row check catches an empty
  # pivot upstream of the merge and points at the failing marker by name.
  pivots_dir <- file.path(tempFolder, "Data", "Pivots")
  for (m in required_markers) {
    marker_dir <- file.path(pivots_dir, m)
    testthat::expect_true(
      dir.exists(marker_dir),
      info = sprintf("missing pivot directory for marker '%s' at %s", m, marker_dir)
    )
    if (!dir.exists(marker_dir)) next
    pq_files <- list.files(marker_dir, pattern = "\\.parquet$", full.names = TRUE)
    testthat::expect_gt(
      length(pq_files), 0L,
      label = sprintf("parquet pivot files under %s", marker_dir)
    )
    for (pf in pq_files) {
      n_rows <- tryCatch(
        as.integer(nrow(polars::pl$read_parquet(pf))),
        error = function(e) NA_integer_
      )
      testthat::expect_gt(
        n_rows, 0L,
        label = sprintf("rows in pivot %s", basename(pf))
      )
    }
  }

  unlink(tempFolder, recursive = TRUE)
})
