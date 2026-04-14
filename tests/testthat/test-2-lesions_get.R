test_that("lesions_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  Sample_ID <- mySampleSheet[1, "Sample_ID"]

  SEMseeker:::init_env(
    result_folder       = tempFolder,
    parallel_strategy   = parallel_strategy,
    maxResources        = 90,
    figures             = "HYPER",
    markers             = "DELTAS",
    areas               = "GENE",
    bonferroni_threshold = 5,
    inpute              = "median"
  )

  if (!exists("signal_thresholds")) {
    signal_data <- SEMseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  values_df <- data.frame(
    CHR   = probe_features$CHR[match(rownames(signal_data), probe_features$PROBE)],
    START = probe_features$START[match(rownames(signal_data), probe_features$PROBE)],
    END   = probe_features$END[match(rownames(signal_data), probe_features$PROBE)],
    VALUE = as.numeric(signal_data[, 1])
  )

  mutations <- SEMseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = signal_thresholds,
    sampleName = Sample_ID
  )

  # ── HYPO lesions: valid data.frame returned ────────────────────────────────
  lesions_hypo <- SEMseeker:::lesions_get(
    mutation_annotated_sorted = mutations,
    grouping_column           = "CHR"
  )

  testthat::expect_s3_class(lesions_hypo, "data.frame")

  # output columns must be exactly CHR / START / END
  testthat::expect_true(all(c("CHR", "START", "END") %in% colnames(lesions_hypo)))

  # lesion count cannot exceed probe count
  testthat::expect_true(nrow(lesions_hypo) <= nrow(mutations))

  # ── HYPER lesions ──────────────────────────────────────────────────────────
  mutations_hyper <- SEMseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = signal_thresholds,
    sampleName = Sample_ID
  )

  lesions_hyper <- SEMseeker:::lesions_get(
    mutation_annotated_sorted = mutations_hyper,
    grouping_column           = "CHR"
  )

  testthat::expect_s3_class(lesions_hyper, "data.frame")

  # ── Edge case: NULL input returns NULL ─────────────────────────────────────
  lesions_null <- SEMseeker:::lesions_get(
    mutation_annotated_sorted = NULL,
    grouping_column           = "CHR"
  )
  testthat::expect_null(lesions_null)

  # ── Edge case: empty data.frame returns 0-row result ──────────────────────
  lesions_empty <- SEMseeker:::lesions_get(
    mutation_annotated_sorted = mutations[0, ],
    grouping_column           = "CHR"
  )
  testthat::expect_equal(nrow(lesions_empty), 0)

  # ── Alternative grouping: PROBE column ────────────────────────────────────
  mutations_with_probe <- mutations
  mutations_with_probe$PROBE <- paste0("cg", formatC(seq_len(nrow(mutations)), width = 7, flag = "0"))
  lesions_probe <- SEMseeker:::lesions_get(
    mutation_annotated_sorted = mutations_with_probe,
    grouping_column           = "PROBE"
  )
  testthat::expect_s3_class(lesions_probe, "data.frame")

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
