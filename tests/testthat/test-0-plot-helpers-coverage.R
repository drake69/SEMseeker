## Coverage for pure helpers behind the plotting layer.
## The plot RENDERERS themselves read pipeline pivots / inference files (or,
## for sem_manhattan_plot_marker_per_sample, reference ssEnv before it is
## assigned), so they need full-pipeline fixtures / a source fix and are not
## exercised here. These two internal helpers are pure and deterministic.
##
## Covered:
##   .assoc_volcano_pick_estimate_col()  choose the estimate column to plot
##   .plot_comparison_pvalue_label()     format a group-comparison p-value label

# ---------------------------------------------------------------------------
# .assoc_volcano_pick_estimate_col
# ---------------------------------------------------------------------------

test_that(".assoc_volcano_pick_estimate_col returns NULL when there is no usable estimate", {
  expect_null(SEMseeker:::.assoc_volcano_pick_estimate_col(data.frame(X = 1, Y = 2)))
  # intercept estimates are excluded
  expect_null(SEMseeker:::.assoc_volcano_pick_estimate_col(
    data.frame(INTERCEPT_ESTIMATE = 1)))
})

test_that(".assoc_volcano_pick_estimate_col falls back to the first estimate w/o PVALUE", {
  expect_equal(
    SEMseeker:::.assoc_volcano_pick_estimate_col(data.frame(AGE_ESTIMATE = 0.5)),
    "AGE_ESTIMATE")
})

test_that(".assoc_volcano_pick_estimate_col matches the estimate whose PVALUE equals PVALUE", {
  df <- data.frame(
    AGE_ESTIMATE = 0.5, AGE_PVALUE = 0.01,
    BMI_ESTIMATE = 0.2, BMI_PVALUE = 0.40,
    PVALUE       = 0.01)
  expect_equal(SEMseeker:::.assoc_volcano_pick_estimate_col(df), "AGE_ESTIMATE")
})

# ---------------------------------------------------------------------------
# .plot_comparison_pvalue_label
# ---------------------------------------------------------------------------

test_that(".plot_comparison_pvalue_label formats a p-value for supported tests", {
  set.seed(1)
  df <- data.frame(
    VALUE = c(stats::rnorm(10, 0), stats::rnorm(10, 4)),
    GRP   = rep(c("a", "b"), each = 10))

  expect_match(
    SEMseeker:::.plot_comparison_pvalue_label(df, "GRP", "VALUE", "t.test"),
    "^p = ")
  expect_match(
    SEMseeker:::.plot_comparison_pvalue_label(df, "GRP", "VALUE", "kruskal.test"),
    "^p = ")
})

test_that(".plot_comparison_pvalue_label returns NA for an unsupported family", {
  df <- data.frame(VALUE = c(1, 2, 3, 4), GRP = c("a", "a", "b", "b"))
  expect_true(is.na(SEMseeker:::.plot_comparison_pvalue_label(df, "GRP", "VALUE", "nope")))
})
