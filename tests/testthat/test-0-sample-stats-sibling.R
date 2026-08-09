## AI-223 slice 1 — per-sample statistics sibling (scope SAMPLE).
##
## Contracts:
##   1. column names come from ONE helper, used by producer and consumer alike;
##   2. the sibling is written next to the sample sheet, one row per sample,
##      carrying burden (SAMPLE_<MARKER>_<FIGURE>) and signal descriptors
##      (SAMPLE_<STAT>);
##   3. the two beta modes are estimated on each side of 0.5 and are omitted on
##      the M-value scale, where the split carries no meaning;
##   4. the sample sheet no longer carries the burden — it moved here;
##   5. the sibling joins back onto the sample sheet on Sample_ID.
##
## Session test uses tempFolders index 13.

# ---------------------------------------------------------------------------
# naming helpers (pure)
# ---------------------------------------------------------------------------

test_that("io_scope_name encodes depth as scope", {
  expect_equal(SEMseeker:::io_scope_name(1L), "SAMPLE")
  expect_equal(SEMseeker:::io_scope_name(1L, "GENE", "BODY"), "SAMPLE")
  expect_equal(SEMseeker:::io_scope_name(2L, "GENE", "BODY"), "GENE")
  expect_equal(SEMseeker:::io_scope_name(3L, "GENE", "BODY"), "GENE_BODY")
  # missing area falls back to sample level whatever the depth
  expect_equal(SEMseeker:::io_scope_name(3L, "", ""), "SAMPLE")
})

test_that("io_stat_colname and io_burden_colname compose the documented names", {
  expect_equal(SEMseeker:::io_stat_colname("SAMPLE", "MEDIAN"), "SAMPLE_MEDIAN")
  expect_equal(SEMseeker:::io_stat_colname("GENE_BODY", "IQR"), "GENE_BODY_IQR")
  expect_equal(
    SEMseeker:::io_burden_colname("SAMPLE", "MUTATIONS", "HYPER"),
    "SAMPLE_MUTATIONS_HYPER")
  expect_equal(
    SEMseeker:::io_burden_colname("GENE_BODY", "LESIONS", "HYPO"),
    "GENE_BODY_LESIONS_HYPO")
  # a scope is mandatory: an unscoped column would be ambiguous in the file
  expect_error(SEMseeker:::io_stat_colname("", "MEDIAN"), "scope")
  expect_error(SEMseeker:::io_burden_colname("", "MUTATIONS", "HYPER"), "scope")
})

test_that("io_signal_stats drops the two modes off the beta scale", {
  expect_true(all(c("MODE_LOW", "MODE_HIGH") %in% SEMseeker:::io_signal_stats(beta = TRUE)))
  expect_false(any(c("MODE_LOW", "MODE_HIGH") %in% SEMseeker:::io_signal_stats(beta = FALSE)))
  expect_true(all(c("MEDIAN", "MEAN", "VARIANCE", "IQR", "N_PROBES") %in%
                    SEMseeker:::io_signal_stats(beta = FALSE)))
})

# ---------------------------------------------------------------------------
# descriptors (pure)
# ---------------------------------------------------------------------------

test_that("the two beta modes land on the injected peaks", {
  set.seed(99L)
  values <- c(stats::rbeta(4000, 2, 40),    # peak near 0.05
              stats::rbeta(4000, 40, 2))    # peak near 0.95

  d <- SEMseeker:::.sem_signal_descriptors(values, beta = TRUE)

  expect_lt(d$MODE_LOW, 0.5)
  expect_gt(d$MODE_HIGH, 0.5)
  expect_lt(d$MODE_LOW, d$MODE_HIGH)
  expect_lt(abs(d$MODE_LOW  - 0.05), 0.1)
  expect_lt(abs(d$MODE_HIGH - 0.95), 0.1)
  expect_equal(d$N_PROBES, length(values))
  expect_equal(d$MEDIAN, stats::median(values))
  expect_equal(d$IQR, stats::IQR(values))
})

test_that("descriptors omit the modes on the M-value scale", {
  set.seed(7L)
  values <- stats::rnorm(1000, mean = 1.5, sd = 2)

  d <- SEMseeker:::.sem_signal_descriptors(values, beta = FALSE)

  expect_null(d$MODE_LOW)
  expect_null(d$MODE_HIGH)
  expect_equal(d$MEAN, mean(values))
  expect_equal(d$VARIANCE, stats::var(values))
})

test_that("descriptors degrade gracefully on degenerate input", {
  empty <- SEMseeker:::.sem_signal_descriptors(numeric(0), beta = TRUE)
  expect_equal(empty$N_PROBES, 0L)
  expect_true(is.na(empty$MEDIAN))

  # a constant vector has no density to speak of
  flat <- SEMseeker:::.sem_signal_descriptors(rep(0.5, 100), beta = TRUE)
  expect_equal(flat$MEDIAN, 0.5)
  expect_true(is.na(flat$MODE_LOW))
  expect_true(is.na(flat$MODE_HIGH))
})

# ---------------------------------------------------------------------------
# end to end
# ---------------------------------------------------------------------------

test_that("semseeker() writes the statistics sibling and leaves the sample sheet lean", {
  tempFolder <- tempFolders[13]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  syn <- .burden_setup_signal_with_outliers()

  SEMseeker::semseeker(
    input             = syn$signal,
    sample_sheet      = syn$samples,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("POSITION"),
    markers           = c("MUTATIONS", "DELTAP", "DELTAS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  stats_csv <- file.path(tempFolder, "Data", "SAMPLE_STATS_RESULT.csv")
  sheet_csv <- file.path(tempFolder, "Data", "SAMPLE_SHEET_RESULT.csv")
  expect_true(file.exists(stats_csv))
  expect_true(file.exists(sheet_csv))
  skip_if_not(file.exists(stats_csv), "downstream assertions need the sibling")

  stats <- utils::read.csv2(stats_csv, stringsAsFactors = FALSE)
  sheet <- utils::read.csv2(sheet_csv, stringsAsFactors = FALSE)

  # one row per sample of the signal matrix
  expect_equal(nrow(stats), ncol(syn$signal))

  descriptor_cols <- c("SAMPLE_MEDIAN", "SAMPLE_MEAN", "SAMPLE_VARIANCE",
                       "SAMPLE_IQR", "SAMPLE_N_PROBES",
                       "SAMPLE_MODE_LOW", "SAMPLE_MODE_HIGH")
  burden_cols <- c("SAMPLE_MUTATIONS_HYPER", "SAMPLE_MUTATIONS_HYPO",
                   "SAMPLE_DELTAP_HYPER", "SAMPLE_DELTAP_HYPO",
                   "SAMPLE_DELTAS_HYPER", "SAMPLE_DELTAS_HYPO")
  expect_true(all(descriptor_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(descriptor_cols, colnames(stats)), collapse = ", ")))
  expect_true(all(burden_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(burden_cols, colnames(stats)), collapse = ", ")))

  # the fixture is on the beta scale: descriptors must be in range and the two
  # modes must sit on either side of 0.5
  expect_true(all(stats$SAMPLE_MEDIAN >= 0 & stats$SAMPLE_MEDIAN <= 1))
  expect_true(all(stats$SAMPLE_IQR >= 0))
  expect_true(all(stats$SAMPLE_N_PROBES > 0))
  expect_true(all(stats$SAMPLE_MODE_LOW < 0.5, na.rm = TRUE))
  expect_true(all(stats$SAMPLE_MODE_HIGH >= 0.5, na.rm = TRUE))

  # AI-223 net move: the burden left the sample sheet
  moved_away <- c("MUTATIONS_HYPER", "MUTATIONS_HYPO", "DELTAS_HYPER",
                  "DELTAP_HYPER", "PROBES_COUNT")
  expect_equal(intersect(moved_away, colnames(sheet)), character(0))

  # and the two files join on Sample_ID
  expect_true(all(sheet$Sample_ID %in% stats$Sample_ID))

  # the join is what consumers see through sem_study_summary_get()
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"),
                            markers = c("MUTATIONS", "DELTAP", "DELTAS"),
                            start_fresh = FALSE, showprogress = FALSE, verbosity = 1)
  joined <- SEMseeker:::sem_study_summary_get()
  expect_true(all(burden_cols %in% colnames(joined)))
  expect_true("SAMPLE_MEDIAN" %in% colnames(joined))
})
