# Regression: position_pivot_to_probe must accept a polars LazyFrame and do
# the join/filter/drop fully in polars (single $collect()), without ever
# materialising the full SIGNAL pivot as an R data.frame. The previous
# implementation R-subset/sorted the collected matrix and produced 4–5 full
# copies of the input, blowing past 80 GB peak on ~367k × 4k inputs and
# triggering macOS Jetsam OOM-kill (Mac/64GB, 2026-06-05).

test_that("position_pivot_to_probe accepts a polars LazyFrame input", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_lazy_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Tiny synthetic annotation: 5 probes on chr1, two samples worth of values.
  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:5)),
    CHR   = "1",
    START = c(100L, 200L, 300L, 400L, 500L),
    END   = c(100L, 200L, 300L, 400L, 500L),
    K850  = TRUE,
    stringsAsFactors = FALSE
  )

  pivot_df <- data.frame(
    CHR   = anno$CHR,
    START = anno$START,
    END   = anno$END,
    S001  = c(0.10, 0.20, 0.30, 0.40, 0.50),
    S002  = c(0.15, 0.25, 0.35, 0.45, 0.55),
    stringsAsFactors = FALSE
  )
  pivot_lazy <- polars::as_polars_df(pivot_df)$lazy()
  testthat::expect_s3_class(pivot_lazy, "polars_lazy_frame")

  # Stub probe_annotation_build so the test does not require the Bioconductor
  # annotation package (macOS setup.R skips loading it).
  testthat::local_mocked_bindings(
    probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  out <- SEMseeker:::position_pivot_to_probe(pivot_lazy)

  testthat::expect_s3_class(out, "data.frame")
  # Output rows are probe-keyed via rownames; no leftover annotation columns.
  testthat::expect_setequal(rownames(out), anno$PROBE)
  testthat::expect_false("PROBE" %in% colnames(out))
  testthat::expect_false("CHR"   %in% colnames(out))
  testthat::expect_false("START" %in% colnames(out))
  testthat::expect_false("END"   %in% colnames(out))
  testthat::expect_false("K850"  %in% colnames(out))
  testthat::expect_setequal(colnames(out), c("S001", "S002"))
  # Values preserved (compare in rowname order to be order-independent).
  testthat::expect_equal(out[anno$PROBE, "S001"], pivot_df$S001)
  testthat::expect_equal(out[anno$PROBE, "S002"], pivot_df$S002)
})

test_that("position_pivot_to_probe still accepts an R data.frame input", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_df_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:3)),
    CHR   = "1",
    START = c(100L, 200L, 300L),
    END   = c(100L, 200L, 300L),
    K850  = TRUE,
    stringsAsFactors = FALSE
  )
  pivot_df <- data.frame(
    CHR = anno$CHR, START = anno$START, END = anno$END,
    S001 = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  out <- SEMseeker:::position_pivot_to_probe(pivot_df)
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_setequal(rownames(out), anno$PROBE)
  testthat::expect_setequal(colnames(out), "S001")
})

test_that("position_pivot_to_probe drops probes flagged FALSE for the tech", {

  skip_on_cran()

  tempFolder <- tempfile("ppp_filter_")
  dir.create(file.path(tempFolder, "Data"), recursive = TRUE)
  ssEnv <- SEMseeker:::init_env(tempFolder,
                                parallel_strategy = "sequential",
                                tech = "K850",
                                iqrTimes = 3, verbosity = 1)
  on.exit({ SEMseeker:::close_env(); unlink(tempFolder, recursive = TRUE) },
          add = TRUE)

  # Two probes ok for K850, one NOT on K850 (should be dropped).
  anno <- data.frame(
    PROBE = paste0("cg", sprintf("%08d", 1:3)),
    CHR   = "1",
    START = c(100L, 200L, 300L),
    END   = c(100L, 200L, 300L),
    K850  = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  pivot_df <- data.frame(
    CHR = anno$CHR, START = anno$START, END = anno$END,
    S001 = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    probe_annotation_build = function(tech) anno,
    .package = "SEMseeker"
  )

  out <- SEMseeker:::position_pivot_to_probe(polars::as_polars_df(pivot_df)$lazy())
  testthat::expect_setequal(rownames(out), anno$PROBE[anno$K850])
  testthat::expect_equal(nrow(out), 2L)
})
