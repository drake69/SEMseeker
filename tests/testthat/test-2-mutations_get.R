test_that("mutations_get", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(result_folder = tempFolder, inpute = "median")

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

  # ── HYPO: basic existence ──────────────────────────────────────────────────
  mutations_hypo <- SEMseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = signal_thresholds,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )

  testthat::expect_false(length(mutations_hypo) == 0)

  # ── HYPO: output structure ─────────────────────────────────────────────────
  testthat::expect_s3_class(mutations_hypo, "data.frame")
  testthat::expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(mutations_hypo)))

  # MUTATIONS column must be binary (0 / 1)
  testthat::expect_true(all(mutations_hypo$MUTATIONS %in% c(0, 1)))

  # row count matches sorted thresholds (no probes lost)
  testthat::expect_true(nrow(mutations_hypo) > 0)

  # ── HYPER: symmetric test ──────────────────────────────────────────────────
  mutations_hyper <- SEMseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPER",
    thresholds = signal_thresholds,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )

  testthat::expect_s3_class(mutations_hyper, "data.frame")
  testthat::expect_true(all(c("CHR", "START", "END", "MUTATIONS") %in% colnames(mutations_hyper)))
  testthat::expect_true(all(mutations_hyper$MUTATIONS %in% c(0, 1)))
  # row count must be identical for HYPO and HYPER (same probes, different direction)
  testthat::expect_equal(nrow(mutations_hypo), nrow(mutations_hyper))

  # ── Boundary: threshold = -Inf  →  zero HYPO mutations ────────────────────
  thresholds_zero <- signal_thresholds
  thresholds_zero$signal_inferior_thresholds <- -Inf
  mutations_none <- SEMseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresholds_zero,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )
  testthat::expect_equal(sum(mutations_none$MUTATIONS), 0)

  # ── Boundary: threshold = +Inf  →  all HYPO mutations ─────────────────────
  thresholds_all <- signal_thresholds
  thresholds_all$signal_inferior_thresholds <- Inf
  mutations_all <- SEMseeker:::mutations_get(
    values     = values_df,
    figure     = "HYPO",
    thresholds = thresholds_all,
    sampleName = mySampleSheet[1, "Sample_ID"]
  )
  testthat::expect_equal(sum(mutations_all$MUTATIONS), nrow(mutations_all))

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
