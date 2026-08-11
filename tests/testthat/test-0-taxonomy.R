## AI-248 — taxonomy SCOPE x MARKER x FIGURE x AGGREGATION.
##
## Contracts:
##   1. exactly two compositors — one for column names (io_feature_colname),
##      one for file names (io_pivot_file_name*); nothing else pastes names;
##   2. the aggregation vocabulary is declared in one place and computed in one
##      place, so a column called MEDIAN cannot hold a mean;
##   3. the FIGURE of SIGNAL is the scale, not an aggregation.

# ---------------------------------------------------------------------------
# column names
# ---------------------------------------------------------------------------

test_that("io_feature_colname carries the four axes in order", {
  expect_equal(
    SEMseeker:::io_feature_colname("SAMPLE", "SIGNAL", "BETA", "MEDIAN"),
    "SAMPLE_SIGNAL_BETA_MEDIAN")
  expect_equal(
    SEMseeker:::io_feature_colname("GENE_TSS1500", "MUTATIONS", "HYPER", "SUM"),
    "GENE_TSS1500_MUTATIONS_HYPER_SUM")
  expect_equal(
    SEMseeker:::io_feature_colname("SAMPLE", "SIGNAL", "MVALUE", "IQR"),
    "SAMPLE_SIGNAL_MVALUE_IQR")
})

test_that("io_feature_colname drops the axes that do not apply", {
  # the number of positions is a property of the scope, not an aggregation of a
  # marker: no marker, no figure
  expect_equal(
    SEMseeker:::io_feature_colname("SAMPLE", aggregation = "N_PROBES"),
    "SAMPLE_N_PROBES")
  expect_equal(
    SEMseeker:::io_feature_colname("GENE_TSS1500", aggregation = "N_PROBES"),
    "GENE_TSS1500_N_PROBES")
  # a scope is mandatory: an unscoped column would be ambiguous in the file
  expect_error(SEMseeker:::io_feature_colname(""), "scope")
  expect_error(SEMseeker:::io_feature_colname(NULL), "scope")
})

# ---------------------------------------------------------------------------
# the aggregation vocabulary
# ---------------------------------------------------------------------------

test_that("the two modes are admissible only on the bounded scale", {
  admissible_beta <- SEMseeker:::util_aggregations_allowed("SIGNAL", "BETA",
                                                           default = FALSE)
  admissible_mval <- SEMseeker:::util_aggregations_allowed("SIGNAL", "MVALUE",
                                                           default = FALSE)
  expect_true(all(c("MODE_LOW", "MODE_HIGH") %in% admissible_beta))
  expect_false(any(c("MODE_LOW", "MODE_HIGH") %in% admissible_mval))
  # an unbounded marker has no bimodal split to speak of either
  expect_false(any(c("MODE_LOW", "MODE_HIGH") %in%
                     SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER",
                                                           default = FALSE)))
})

test_that("a 0/1 marker admits only the two aggregations that carry information", {
  admissible <- SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER",
                                                      discrete = TRUE, default = FALSE)
  # the burden, and the burden per position — the density, which is what makes
  # regions of different size comparable
  expect_setequal(admissible, c("SUM", "MEAN"))
  # median / variance / IQR of a vector of zeros and ones are degenerate
  expect_false(any(c("MEDIAN", "VARIANCE", "IQR") %in% admissible))
})

test_that("the continuous markers are one class: DELTAS and DELTAR match SIGNAL", {
  signal <- SEMseeker:::util_aggregations_allowed("SIGNAL", "MVALUE",
                                                  discrete = FALSE, default = FALSE)
  deltas <- SEMseeker:::util_aggregations_allowed("DELTAS", "HYPER",
                                                  discrete = FALSE, default = FALSE)
  deltar <- SEMseeker:::util_aggregations_allowed("DELTAR", "HYPO",
                                                  discrete = FALSE, default = FALSE)
  expect_setequal(deltas, signal)
  expect_setequal(deltar, signal)
  expect_true(all(c("MEAN", "MEDIAN", "VARIANCE", "IQR") %in% deltas))
})

test_that("what is produced by default is narrower than what is admissible", {
  # a count of epimutations is a burden, and the burden is its sum
  expect_equal(SEMseeker:::util_aggregations_allowed("MUTATIONS", "HYPER"), "SUM")
  # every continuous marker carries the whole descriptor set: adding them does
  # not change the historical mean, it only gives it a suffix
  for (marker in c("SIGNAL", "DELTAS", "DELTAR")) {
    produced <- SEMseeker:::util_aggregations_allowed(marker, "HYPER",
                                                      discrete = FALSE)
    expect_true(all(c("MEAN", "MEDIAN", "VARIANCE", "IQR") %in% produced),
                info = marker)
    expect_false("SUM" %in% produced, info = marker)
  }
})

# ---------------------------------------------------------------------------
# what each aggregation computes
# ---------------------------------------------------------------------------

test_that("util_aggregate_values computes what the name says", {
  set.seed(3L)
  values <- stats::rnorm(500, mean = 2, sd = 3)

  expect_equal(SEMseeker:::util_aggregate_values(values, "SUM"), sum(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "MEAN"), mean(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "MEDIAN"),
               stats::median(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "VARIANCE"),
               stats::var(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "IQR"), stats::IQR(values))
  expect_equal(SEMseeker:::util_aggregate_values(values, "N_PROBES"), length(values))
})

test_that("util_aggregate_values finds the two peaks of a bimodal sample", {
  set.seed(11L)
  values <- c(stats::rbeta(3000, 2, 40), stats::rbeta(3000, 40, 2))
  low  <- SEMseeker:::util_aggregate_values(values, "MODE_LOW")
  high <- SEMseeker:::util_aggregate_values(values, "MODE_HIGH")
  expect_lt(low, 0.5)
  expect_gt(high, 0.5)
})

test_that("util_aggregate_values degrades instead of inventing a number", {
  expect_true(is.na(SEMseeker:::util_aggregate_values(numeric(0), "MEDIAN")))
  expect_equal(SEMseeker:::util_aggregate_values(numeric(0), "N_PROBES"), 0L)
  # non-finite values are not data
  expect_equal(SEMseeker:::util_aggregate_values(c(1, NA, Inf, 3), "N_PROBES"), 2L)
  expect_error(SEMseeker:::util_aggregate_values(1:10, "AVERAGE"), "unknown aggregation")
})

# ---------------------------------------------------------------------------
# the figure of SIGNAL is the scale
# ---------------------------------------------------------------------------

test_that("io_signal_figure reports the scale, not an aggregation", {
  expect_equal(SEMseeker:::io_signal_figure(beta = TRUE), "BETA")
  expect_equal(SEMseeker:::io_signal_figure(beta = FALSE), "MVALUE")
  # a session that has not seen the data yet defaults to the bounded scale
  expect_equal(SEMseeker:::io_signal_figure(beta = NA), "BETA")
})

# ---------------------------------------------------------------------------
# file names
# ---------------------------------------------------------------------------

test_that("the pivot name carries the aggregation, and the scale for SIGNAL", {
  tempFolder <- tempFolders[17]
  unlink(tempFolder, recursive = TRUE)
  on.exit({ try(SEMseeker:::core_close_env(), silent = TRUE)
            unlink(tempFolder, recursive = TRUE) }, add = TRUE)

  SEMseeker:::core_init_env(result_folder = tempFolder, parallel_strategy = "sequential",
                            areas = c("POSITION"), markers = c("MUTATIONS"),
                            start_fresh = TRUE, showprogress = FALSE, verbosity = 1)

  burden <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER", "CHR",
                                                   "CYTOBAND", aggregation = "SUM")
  expect_match(basename(burden), "^MUTATIONS_HYPER_CHR_CYTOBAND_SUM_HG19\\.parquet$",
               ignore.case = TRUE)

  beta_mean <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "BETA", "CHR",
                                                      "CYTOBAND", aggregation = "MEAN")
  mval_mean <- SEMseeker:::io_pivot_file_name_parquet("SIGNAL", "MVALUE", "CHR",
                                                      "CYTOBAND", aggregation = "MEAN")
  # the two scales must not overwrite each other in the same folder
  expect_false(identical(beta_mean, mval_mean))

  # omitting the aggregation is how the not-yet-aggregated POSITION pivots are
  # named, and must keep the historical shape
  position <- SEMseeker:::io_pivot_file_name_parquet("MUTATIONS", "HYPER",
                                                     "POSITION", "WHOLE")
  expect_match(basename(position), "^MUTATIONS_HYPER_POSITION_WHOLE_HG19\\.parquet$",
               ignore.case = TRUE)
})

# ---------------------------------------------------------------------------
# coherence of the request, checked at the door
# ---------------------------------------------------------------------------

.tax_keys <- function(...) {
  k <- data.frame(..., stringsAsFactors = FALSE)
  k
}

test_that("a request that names no aggregation is refused, at every depth", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  for (depth in c(1L, 2L, 3L)) {
    details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                          depth_analysis = depth, aggregation = NA,
                          stringsAsFactors = FALSE)
    expect_error(SEMseeker:::assoc_validate_aggregation(details, keys), "required")
  }
})

test_that("an aggregation outside the taxonomy is refused", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                        depth_analysis = 1L, aggregation = "AVERAGE",
                        stringsAsFactors = FALSE)
  expect_error(SEMseeker:::assoc_validate_aggregation(details, keys),
               "not an aggregation")
})

test_that("an impossible request drops its row instead of stopping the batch", {
  keys <- .tax_keys(MARKER = c("MUTATIONS", "SIGNAL"),
                    FIGURE = c("HYPER", "BETA"),
                    DISCRETE = c(TRUE, FALSE))
  details <- data.frame(
    independent_variable = c("Phenotest", "Phenotest"),
    family_test          = c("spearman", "spearman"),
    depth_analysis       = c(1L, 1L),
    # MODE_LOW exists only for SIGNAL/BETA, so the first row is possible;
    # a run with only counts would make it impossible
    aggregation          = c("MEDIAN", "SUM"),
    stringsAsFactors = FALSE)

  expect_warning(kept <- SEMseeker:::assoc_validate_aggregation(details, keys),
                 "not admissible for MUTATIONS/HYPER")
  # partially admissible: the row survives, the impossible pairs are skipped
  expect_equal(nrow(kept), 2L)

  # now the same median against counts only: nothing can be computed
  counts_only <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  one_row <- details[1, , drop = FALSE]
  expect_error(
    suppressWarnings(SEMseeker:::assoc_validate_aggregation(one_row, counts_only)),
    "nothing left to analyse")
})

test_that("a batch keeps the rows it can run and drops the ones it cannot", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  details <- data.frame(
    independent_variable = c("Phenotest", "Phenotest"),
    family_test          = c("spearman", "spearman"),
    depth_analysis       = c(1L, 1L),
    aggregation          = c("MEDIAN", "SUM"),
    stringsAsFactors = FALSE)

  expect_warning(kept <- SEMseeker:::assoc_validate_aggregation(details, keys),
                 "row 1")
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$aggregation, "SUM")
})

test_that("a possible request passes silently", {
  keys <- .tax_keys(MARKER = "MUTATIONS", FIGURE = "HYPER", DISCRETE = TRUE)
  details <- data.frame(independent_variable = "Phenotest", family_test = "spearman",
                        depth_analysis = 1L, aggregation = "SUM",
                        stringsAsFactors = FALSE)
  expect_silent(kept <- SEMseeker:::assoc_validate_aggregation(details, keys))
  expect_equal(nrow(kept), 1L)
})
