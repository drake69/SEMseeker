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

test_that("io_scope_name derives the scope from (area, subarea) with no depth", {
  # AI-223 slice 2a: the producer is depth-agnostic — it writes every scope
  # once, so it names them from the pair alone.
  expect_equal(SEMseeker:::io_scope_name(area = "GENE", subarea = "TSS1500"),
               "GENE_TSS1500")
  expect_equal(SEMseeker:::io_scope_name(area = "GENE"), "GENE")
  expect_equal(SEMseeker:::io_scope_name(), "SAMPLE")
  expect_equal(SEMseeker:::io_scope_name(area = "", subarea = "TSS1500"), "SAMPLE")
})

# ---------------------------------------------------------------------------
# descriptors (pure)
# ---------------------------------------------------------------------------

test_that("the two beta modes land on the injected peaks", {
  set.seed(99L)
  values <- c(stats::rbeta(4000, 2, 40),    # peak near 0.05
              stats::rbeta(4000, 40, 2))    # peak near 0.95

  d <- SEMseeker:::util_signal_descriptors(values, beta = TRUE)

  expect_lt(d$MODELOW, 0.5)
  expect_gt(d$MODEHIGH, 0.5)
  expect_lt(d$MODELOW, d$MODEHIGH)
  expect_lt(abs(d$MODELOW  - 0.05), 0.1)
  expect_lt(abs(d$MODEHIGH - 0.95), 0.1)
  expect_equal(d$N_PROBES, length(values))
  expect_equal(d$MEDIAN, stats::median(values))
  expect_equal(d$IQR, stats::IQR(values))
})

test_that("descriptors omit the modes on the M-value scale", {
  set.seed(7L)
  values <- stats::rnorm(1000, mean = 1.5, sd = 2)

  d <- SEMseeker:::util_signal_descriptors(values, beta = FALSE)

  expect_null(d$MODELOW)
  expect_null(d$MODEHIGH)
  expect_equal(d$MEAN, mean(values))
  expect_equal(d$VARIANCE, stats::var(values))
})

test_that("descriptors degrade gracefully on degenerate input", {
  empty <- SEMseeker:::util_signal_descriptors(numeric(0), beta = TRUE)
  expect_equal(empty$N_PROBES, 0L)
  expect_true(is.na(empty$MEDIAN))

  # a constant vector has no density to speak of
  flat <- SEMseeker:::util_signal_descriptors(rep(0.5, 100), beta = TRUE)
  expect_equal(flat$MEDIAN, 0.5)
  expect_true(is.na(flat$MODELOW))
  expect_true(is.na(flat$MODEHIGH))
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

  # AI-255: the statistics are SCOPE = SAMPLE artefacts, composed on read.
  sheet_csv <- file.path(tempFolder, "Data", "SAMPLE_SHEET_RESULT.csv")
  expect_true(file.exists(sheet_csv))

  stats <- SEMseeker:::sem_study_summary_get()
  sheet <- utils::read.csv2(sheet_csv, stringsAsFactors = FALSE)
  skip_if_not(!is.null(stats) && nrow(stats) > 0,
              "downstream assertions need the composed statistics")

  # one row per sample of the signal matrix
  expect_equal(nrow(stats), ncol(syn$signal))

  # AI-248: `markers` means one thing for every marker, SIGNAL included — this
  # run did not ask for it, so its descriptors must NOT be there. N_PROBES is
  # not a marker column: it is the size of the scope, the denominator of the
  # density, and it is there regardless.
  descriptor_cols <- "SAMPLE_N_PROBES"
  signal_cols <- vapply(c("MEDIAN", "MEAN", "VARIANCE", "IQR", "MODELOW", "MODEHIGH"),
                        function(a) SEMseeker:::io_feature_colname("SAMPLE", "SIGNAL", "BETA", a),
                        character(1))
  burden_cols <- c(
    SEMseeker:::io_feature_colname("SAMPLE", "MUTATIONS", c("HYPER", "HYPO"), "SUM"),
    SEMseeker:::io_feature_colname("SAMPLE", "DELTAP", c("HYPER", "HYPO"), "SUM"),
    SEMseeker:::io_feature_colname("SAMPLE", "DELTAS", c("HYPER", "HYPO"), "MEAN"))
  expect_true(all(descriptor_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(descriptor_cols, colnames(stats)), collapse = ", ")))
  expect_true(all(burden_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(burden_cols, colnames(stats)), collapse = ", ")))

  # a marker that was not asked for produces no column
  expect_equal(intersect(signal_cols, colnames(stats)), character(0))
  expect_true(all(stats$SAMPLE_N_PROBES > 0))

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
  expect_true("SAMPLE_N_PROBES" %in% colnames(joined))
})

# ---------------------------------------------------------------------------
# AI-223 slice 2a — region scope, produced and consumed at depth = 1
# ---------------------------------------------------------------------------

test_that("a region scope reaches the sibling and the depth=1 inference", {
  tempFolder <- tempFolders[14]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  syn <- .burden_setup_signal_with_outliers()
  scope <- SEMseeker:::io_scope_name(area = "GENE", subarea = "TSS1500")

  SEMseeker::semseeker(
    input               = syn$signal,
    sample_sheet        = syn$samples,
    result_folder       = tempFolder,
    parallel_strategy   = "sequential",
    # "WHOLE" must stay in `subareas`: it is what keeps POSITION_WHOLE in the
    # keys, and every burden — whatever its scope — is aggregated from the
    # POSITION pivots.
    areas               = c("POSITION", "GENE"),
    subareas            = c("WHOLE", "TSS1500"),
    markers             = c("MUTATIONS", "SIGNAL"),
    start_fresh         = TRUE,
    inpute              = "median",
    showprogress        = showprogress,
    verbosity           = verbosity
  )

  # AI-255: the region class is asked for at read time, not declared before the
  # run. This is the whole point — no rerun to change your mind.
  stats <- SEMseeker:::sem_study_summary_get(regions = c("SAMPLE", scope))
  skip_if_not(!is.null(stats) && nrow(stats) > 0,
              "downstream assertions need the composed statistics")

  scope_cols <- SEMseeker:::io_feature_colname(scope, "MUTATIONS", c("HYPER", "HYPO"), "SUM")
  whole_cols <- SEMseeker:::io_feature_colname("SAMPLE", "MUTATIONS", c("HYPER", "HYPO"), "SUM")
  expect_true(all(scope_cols %in% colnames(stats)),
              info = paste("missing:",
                           paste(setdiff(scope_cols, colnames(stats)), collapse = ", ")))
  expect_true(all(whole_cols %in% colnames(stats)))

  # a subset of the probes can only carry a subset of the burden
  for (i in seq_along(scope_cols))
    expect_true(all(stats[[scope_cols[i]]] <= stats[[whole_cols[i]]]))

  # AI-248: asking for SIGNAL gives its descriptors on the scope too — this is
  # the per-sample median restricted to a region class
  scope_median <- SEMseeker:::io_feature_colname(scope, "SIGNAL", "BETA", "MEDIAN")
  expect_true(scope_median %in% colnames(stats),
              info = paste("missing:", scope_median))
  expect_true(all(stats[[scope_median]] >= 0 & stats[[scope_median]] <= 1, na.rm = TRUE))
  # the scope holds fewer positions than the whole sample
  expect_true(all(stats[[SEMseeker:::io_feature_colname(scope, aggregation = "N_PROBES")]] <=
                    stats$SAMPLE_N_PROBES))

  # ── the mask counts every position once ──────────────────────────────────
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION", "GENE"), subareas = c("WHOLE", "TSS1500"),
                            markers = c("MUTATIONS"),
                            start_fresh = FALSE, showprogress = FALSE, verbosity = 1)

  mask <- SEMseeker:::.sem_scope_probe_mask("GENE", "TSS1500")
  expect_false(is.null(mask))
  skip_if(is.null(mask), "no probe of the fixture is annotated TSS1500")
  mask_df <- as.data.frame(mask$collect())
  expect_gt(nrow(mask_df), 0)

  pivot <- SEMseeker:::io_read_pivot("MUTATIONS", "HYPER", "POSITION", "WHOLE")
  pivot_df <- as.data.frame(pivot$collect())
  pivot_df$CHR   <- sub("^(?i)chr", "", as.character(pivot_df$CHR), perl = TRUE)
  pivot_df$START <- as.integer(pivot_df$START)
  pivot_df$END   <- as.integer(pivot_df$END)
  masked <- merge(pivot_df, mask_df, by = c("CHR", "START", "END"))
  sample_columns <- setdiff(colnames(masked), c("CHR", "START", "END"))
  expected <- colSums(masked[, sample_columns, drop = FALSE], na.rm = TRUE)

  observed <- stats[[SEMseeker:::io_feature_colname(scope, "MUTATIONS", "HYPER", "SUM")]]
  names(observed) <- stats$Sample_ID
  common <- intersect(names(expected), names(observed))
  expect_gt(length(common), 0)
  expect_equal(as.numeric(observed[common]), as.numeric(expected[common]))

  # a probe annotated to several genes must not be counted twice: the mask has
  # one row per position, whatever the annotation says
  expect_equal(nrow(mask_df), nrow(unique(mask_df[, c("CHR", "START", "END")])))

  # ── consumption at depth = 1 ─────────────────────────────────────────────
  SEMseeker:::core_close_env()

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    depth_analysis       = 1L,
    scopes               = paste("SAMPLE", scope, sep = "+"),
    aggregation          = "SUM",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  expect_no_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      # depth=1 reads the sibling, not the area pivots: no need to pay for the
      # GENE keys here.
      areas             = c("POSITION"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    )
  )

  csv_files <- list.files(file.path(tempFolder, "Inference"), pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  csv_files <- csv_files[!grepl("(?i)assoc_covariates_model", csv_files)]
  expect_gt(length(csv_files), 0)
  skip_if(length(csv_files) == 0, "no inference CSV to inspect")

  result_df <- do.call(plyr::rbind.fill,
                       lapply(csv_files, function(f) utils::read.csv2(f, stringsAsFactors = FALSE)))
  # the scope is one value per sample: still DEPTH 1, with the scope in AREA.
  # which() and not a bare logical: the CSV carries rows with NA in AREA (the
  # job-summary row), and `df[df$AREA == scope, ]` would materialise them as
  # all-NA phantom rows.
  scope_rows <- result_df[which(result_df$AREA == scope), , drop = FALSE]
  expect_gt(nrow(scope_rows), 0)
  expect_true(all(scope_rows$DEPTH == 1))
  expect_true(all(scope_rows$SUBAREA == "SAMPLE"))
  # which burden was tested stays in MARKER/FIGURE, as at every other depth;
  # the scope is what AREA adds. (AREA_OF_TEST is "burden_values" on every
  # depth=1 row, scoped or not: with a single tested column io_data_preparation
  # loses the column name — pre-existing, see AI-247.)
  expect_equal(sort(unique(scope_rows$MARKER)), "MUTATIONS")
  expect_true(all(scope_rows$FIGURE %in% c("HYPER", "HYPO")))
  expect_setequal(unique(scope_rows$FIGURE),
                  unique(result_df[which(result_df$AREA == "SAMPLE_GROUP"), "FIGURE"]))
  # and the whole-sample rows are still there, untouched
  expect_gt(nrow(result_df[which(result_df$AREA == "SAMPLE_GROUP"), , drop = FALSE]), 0)
})

test_that("an unproduced scope stops the analysis instead of testing nothing", {
  tempFolder <- tempFolders[15]
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
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    inpute            = "median",
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    depth_analysis       = 1L,
    scopes               = "GENE_TSS1500",
    aggregation          = "SUM",
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  expect_error(
    SEMseeker:::association_analysis(
      inference_details = inference_details,
      result_folder     = tempFolder,
      parallel_strategy = "sequential",
      markers           = c("MUTATIONS"),
      areas             = c("POSITION"),
      multiple_test_adj = "BH",
      showprogress      = showprogress,
      verbosity         = verbosity
    ),
    "GENE_TSS1500")
})

test_that("an unknown region class is refused at the door, not silently ignored", {
  tempFolder <- tempFolders[16]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  # AI-255: the region classes are no longer declared at SEM time — semseeker()
  # lost sample_stats_scopes — so the refusal moved to where they are asked for.
  # A run that quietly dropped an unknown class would look identical to one that
  # honoured it, and the researcher would find out only at analysis time.
  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"), markers = c("MUTATIONS"),
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  expect_error(SEMseeker:::.sem_regions_resolve(c("SAMPLE", "GENE_NOWHERE")),
               "GENE_NOWHERE")
  # ... and the one that is not a region class at all still resolves, because
  # "no restriction" is a legal request.
  expect_equal(SEMseeker:::.sem_regions_resolve("SAMPLE")[[1]]$area, "PROBE")
})
