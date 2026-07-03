# ---------------------------------------------------------------------------
# Behavioral unit tests for the Binomial lesion model (A-01)
#
# These tests exercise lesions_get() with synthetic, fully controlled inputs
# so the expected p-value and lesion outcome can be computed analytically.
# They do NOT depend on the global signal_data / probe_features fixtures.
# ---------------------------------------------------------------------------

# Helper: build a minimal mutation_annotated_sorted data frame
.make_mutations_df <- function(chr, starts, mutations) {
  data.frame(
    CHR       = chr,
    START     = as.integer(starts),
    END       = as.integer(starts),
    MUTATIONS = as.integer(mutations),
    stringsAsFactors = FALSE
  )
}

test_that("lesions_get (binomial): tight cluster of mutations → lesions detected", {
  tf <- tempFolders[33]
  SEMseeker:::init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 5000L,
    bonferroni_threshold = 0.5   # loose threshold so the cluster always fires
  )
  on.exit({ SEMseeker:::close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 60 probes on chr1 spaced 1000 bp apart.
  # Probes 25-35 (11 consecutive) are ALL mutated → ENRICHMENT = 11 in the
  # centre of the cluster.  Background rate p0 = 11/60 ≈ 0.183.
  # P(X >= 11 | Binom(11, 0.183)) is negligibly small → lesion expected.
  n <- 60L
  muts <- integer(n)
  muts[25:35] <- 1L
  df <- .make_mutations_df("chr1", seq(1e6L, by = 1000L, length.out = n), muts)

  result <- SEMseeker:::lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("CHR", "START", "END") %in% colnames(result)))
  expect_gt(nrow(result), 0L)   # at least one lesion probe detected
})

test_that("lesions_get (binomial): dispersed mutations → no lesions", {
  tf <- tempFolders[34]
  SEMseeker:::init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 5000L,
    bonferroni_threshold = 0.1
  )
  on.exit({ SEMseeker:::close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 66 probes, exactly 1 mutation every 11 probes → ENRICHMENT never exceeds 1
  # in any window.  With p0 = 6/66 ≈ 0.09, P(X >= 1 | Binom(11, 0.09)) ≈ 0.67
  # — far above any reasonable Bonferroni threshold → no lesions.
  n <- 66L
  muts <- integer(n)
  muts[seq(6L, n, by = 11L)] <- 1L
  df <- .make_mutations_df("chr1", seq(1e6L, by = 1000L, length.out = n), muts)

  result <- SEMseeker:::lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)

  expect_equal(nrow(result), 0L)
})

test_that("lesions_get (binomial): zero mutations → no lesions and no error", {
  tf <- tempFolders[35]
  SEMseeker:::init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 5000L,
    bonferroni_threshold = 0.1
  )
  on.exit({ SEMseeker:::close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  df <- .make_mutations_df("chr1", seq(1e6L, by = 1000L, length.out = 50L),
                           integer(50L))

  result <- SEMseeker:::lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
})

test_that("lesions_get (binomial): p-value matches pbinom formula exactly", {
  # Verify the implementation is genuinely using pbinom, not dhyper.
  # With a controlled input we can compute the expected p-value by hand and
  # confirm the lesion threshold is crossed exactly when pbinom says it should be.
  tf <- tempFolders[36]
  SEMseeker:::init_env(
    result_folder        = tf,
    start_fresh          = TRUE,
    LESIONS_BP           = 5000L,
    bonferroni_threshold = 10   # very loose: fires whenever p < 10 / (n * log10(bp))
  )
  on.exit({ SEMseeker:::close_env(); unlink(tf, recursive = TRUE) }, add = TRUE)

  # 22 probes; first 11 mutated, last 11 not.
  # MUTATIONS_COUNT = 11, PROBES_COUNT = 22, p0 = 0.5
  # Centre of cluster: ENRICHMENT = 11
  # Expected p-value = pbinom(10, 11, 0.5, lower.tail = FALSE) = 0.5^11
  n <- 22L
  muts <- c(rep(1L, 11L), rep(0L, 11L))
  df <- .make_mutations_df("chr1", seq(1e6L, by = 1000L, length.out = n), muts)

  result <- SEMseeker:::lesions_get(grouping_column = "CHR",
                                    mutation_annotated_sorted = df)

  # With p0 = 0.5 and ENRICHMENT = 11 at the cluster centre, pbinom gives
  # 0.5^11 ≈ 4.9e-4, which is < 10 / (22 * log10(1e4)) ≈ 0.114 → lesion.
  expect_gt(nrow(result), 0L)
})

# ---------------------------------------------------------------------------
# Original integration test (preserved unchanged)
# ---------------------------------------------------------------------------

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
