# Tests for pure (side-effect-free) helper functions
#
# Covered:
#  - name_cleaning()                      string sanitiser
#  - is.family_dicotomic()                family-test classifier
#  - sig.formula_vars()                   formula parser
#  - apply_stat_model_sig_formula()       formula builder
#  - exact_pvalue()                       permutation p-value calculator
#  - model_performance()                  regression metric calculator
#  - sort_by_chr_and_start()              genomic-coordinate sorter
#  - compute_mean_delta_permutation_cpu() mean-difference permutation statistic
#  - compute_spearman_permutation()       Spearman permutation statistic
#  - compute_quantile_delta_permutation() quantile-difference permutation statistic
#  - describe_dataframe()                 per-column summary statistics
#  - polynomial_formula_build()           polynomial regression formula builder

# ---------------------------------------------------------------------------
# 1. name_cleaning
# ---------------------------------------------------------------------------

test_that("name_cleaning uppercases and replaces spaces/dots/dashes with underscores", {
  expect_equal(SEMseeker:::name_cleaning("hello world"), "HELLO_WORLD")
  expect_equal(SEMseeker:::name_cleaning("a.b.c"),       "A_B_C")
  expect_equal(SEMseeker:::name_cleaning("a-b-c"),       "A_B_C")
  expect_equal(SEMseeker:::name_cleaning("abc"),         "ABC")
})

test_that("name_cleaning replaces comparison operators", {
  expect_equal(SEMseeker:::name_cleaning("x>=1"),  "X_GTE_1")
  expect_equal(SEMseeker:::name_cleaning("x<=1"),  "X_LTE_1")
  expect_equal(SEMseeker:::name_cleaning("x==1"),  "X_EQ_1")
  expect_equal(SEMseeker:::name_cleaning("x!=1"),  "X_DISEQ_1")
  expect_equal(SEMseeker:::name_cleaning("x>1"),   "X_GT_1")
  expect_equal(SEMseeker:::name_cleaning("x<1"),   "X_LT_1")
})

test_that("name_cleaning strips leading and trailing underscores", {
  expect_false(startsWith(SEMseeker:::name_cleaning("_abc"), "_"))
  expect_false(endsWith(SEMseeker:::name_cleaning("abc_"),   "_"))
})

test_that("name_cleaning collapses double underscores", {
  result <- SEMseeker:::name_cleaning("a__b")
  expect_false(grepl("__", result))
})

test_that("name_cleaning handles empty strings with placeholder", {
  expect_equal(SEMseeker:::name_cleaning("", empty_place_holder = "EMPTY"), "EMPTY")
})

test_that("name_cleaning handles vectors", {
  result <- SEMseeker:::name_cleaning(c("hello world", "foo.bar"))
  expect_equal(result, c("HELLO_WORLD", "FOO_BAR"))
})

# ---------------------------------------------------------------------------
# 2. is.family_dicotomic
# ---------------------------------------------------------------------------

test_that("is.family_dicotomic returns FALSE for continuous families", {
  expect_false(SEMseeker:::is.family_dicotomic("gaussian"))
  expect_false(SEMseeker:::is.family_dicotomic("spearman"))
  expect_false(SEMseeker:::is.family_dicotomic("pearson"))
  expect_false(SEMseeker:::is.family_dicotomic("kendall"))
  expect_false(SEMseeker:::is.family_dicotomic("poisson"))
  expect_false(SEMseeker:::is.family_dicotomic("quantreg_0.5"))
  expect_false(SEMseeker:::is.family_dicotomic("quantreg-permutation_0.5_5_10_0.9"))
  expect_false(SEMseeker:::is.family_dicotomic("mean-permutation_100_1000_0.95"))
  expect_false(SEMseeker:::is.family_dicotomic("spearman-permutation"))
  expect_false(SEMseeker:::is.family_dicotomic("quantile-permutation"))
  expect_false(SEMseeker:::is.family_dicotomic("polynomial_4_1"))
  expect_false(SEMseeker:::is.family_dicotomic("log_1"))
})

test_that("is.family_dicotomic returns TRUE for dichotomous families", {
  expect_true(SEMseeker:::is.family_dicotomic("wilcoxon"))
  expect_true(SEMseeker:::is.family_dicotomic("t.test"))
  expect_true(SEMseeker:::is.family_dicotomic("kruskal.test"))
  expect_true(SEMseeker:::is.family_dicotomic("chisq.test"))
  expect_true(SEMseeker:::is.family_dicotomic("fisher.test"))
  expect_true(SEMseeker:::is.family_dicotomic("binomial"))
  expect_true(SEMseeker:::is.family_dicotomic("multinomial"))
  expect_true(SEMseeker:::is.family_dicotomic("jsd"))
})

# ---------------------------------------------------------------------------
# 3. sig.formula_vars
# ---------------------------------------------------------------------------

test_that("sig.formula_vars extracts dependent and independent variables", {
  f <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::sig.formula_vars(f)
  expect_equal(res$dependent_variable,   "burden")
  expect_equal(res$independent_variable, "group")
  expect_length(res$covariates, 0)
})

test_that("sig.formula_vars extracts covariates correctly", {
  f <- stats::as.formula("burden ~ group + age + sex")
  res <- SEMseeker:::sig.formula_vars(f)
  expect_equal(res$dependent_variable,   "burden")
  expect_equal(res$independent_variable, "group")
  expect_equal(sort(res$covariates), c("age", "sex"))
})

# ---------------------------------------------------------------------------
# 4. apply_stat_model_sig_formula
# ---------------------------------------------------------------------------

test_that("apply_stat_model_sig_formula builds correct formula for wilcoxon", {
  f <- SEMseeker:::apply_stat_model_sig_formula("wilcoxon", "BURDEN", "GROUP", NULL)
  expect_equal(deparse(f), "BURDEN ~ GROUP")
})

test_that("apply_stat_model_sig_formula builds correct formula for gaussian with covariates", {
  f <- SEMseeker:::apply_stat_model_sig_formula("gaussian", "BURDEN", "GROUP", c("AGE", "SEX"))
  expect_true(grepl("GROUP", deparse(f)))
  expect_true(grepl("AGE",   deparse(f)))
  expect_true(grepl("SEX",   deparse(f)))
})

test_that("apply_stat_model_sig_formula inverts roles for binomial", {
  f <- SEMseeker:::apply_stat_model_sig_formula("binomial", "BURDEN", "GROUP", NULL)
  # for binomial: GROUP ~ BURDEN
  expect_equal(deparse(f), "GROUP ~ BURDEN")
})

test_that("apply_stat_model_sig_formula builds correct formula for pearson", {
  f <- SEMseeker:::apply_stat_model_sig_formula("pearson", "BURDEN", "PHENO", NULL)
  expect_equal(deparse(f), "BURDEN ~ PHENO")
})

test_that("apply_stat_model_sig_formula builds correct formula for mean-permutation", {
  f <- SEMseeker:::apply_stat_model_sig_formula("mean-permutation_100_1000_0.95", "CG_AREA", "SAMPLE_GROUP", NULL)
  expect_equal(deparse(f), "CG_AREA ~ SAMPLE_GROUP")
})

# ---------------------------------------------------------------------------
# 5. exact_pvalue
# ---------------------------------------------------------------------------

test_that("exact_pvalue returns vector of length 4", {
  set.seed(42)
  perm <- rnorm(1000, mean = 0, sd = 1)
  result <- SEMseeker:::exact_pvalue(perm, statistic_parameter = 3, conf.level = 0.95)
  expect_length(result, 4)
})

test_that("exact_pvalue: large statistic gives small p-value", {
  set.seed(42)
  perm <- rnorm(10000, mean = 0, sd = 1)
  result <- SEMseeker:::exact_pvalue(perm, statistic_parameter = 5, conf.level = 0.95)
  p_value <- result[3]
  expect_lt(p_value, 0.01)
})

test_that("exact_pvalue: statistic = median gives p-value ~ 0.5", {
  set.seed(42)
  perm <- rnorm(10000, mean = 0, sd = 1)
  med <- median(perm)
  result <- SEMseeker:::exact_pvalue(perm, statistic_parameter = med, conf.level = 0.95)
  p_value <- result[3]
  expect_lt(p_value, 0.6)
  expect_gt(p_value, 0.4)
})

test_that("exact_pvalue: CI lower < CI upper", {
  set.seed(42)
  perm <- rnorm(1000)
  result <- SEMseeker:::exact_pvalue(perm, statistic_parameter = 0, conf.level = 0.95)
  expect_lt(result[1], result[2])
})

test_that("exact_pvalue: pvalue_limit = 1 - conf.level", {
  set.seed(42)
  perm <- rnorm(1000)
  conf_level <- 0.90
  result <- SEMseeker:::exact_pvalue(perm, statistic_parameter = 0, conf.level = conf_level)
  expect_equal(unname(result[4]), 1 - conf_level)
})

# ---------------------------------------------------------------------------
# 6. model_performance
# ---------------------------------------------------------------------------

test_that("model_performance returns a data.frame with expected columns", {
  fitted   <- c(1.1, 2.0, 3.1, 4.0)
  expected <- c(1.0, 2.0, 3.0, 4.0)
  res <- SEMseeker:::model_performance(fitted, expected, c(), c())
  expect_s3_class(res, "data.frame")
  expect_true(all(c("mse", "rmse", "mae", "r_squared") %in% colnames(res)))
})

test_that("model_performance: perfect fit gives r_squared = 1", {
  x <- c(1, 2, 3, 4, 5)
  res <- SEMseeker:::model_performance(x, x, c(), c())
  expect_equal(res$mse, 0)
  expect_equal(res$mae, 0)
  expect_equal(res$rmse, 0)
  expect_equal(res$r_squared, 1)
})

test_that("model_performance: mse is mean of squared residuals", {
  fitted   <- c(2, 4)
  expected <- c(1, 3)
  res <- SEMseeker:::model_performance(fitted, expected, c(), c())
  expect_equal(res$mse, 1)   # both residuals = -1, mean((1)^2) = 1
})

test_that("model_performance: rmse = sqrt(mse)", {
  fitted   <- c(2, 2, 6)
  expected <- c(1, 3, 5)
  res <- SEMseeker:::model_performance(fitted, expected, c(), c())
  expect_equal(res$rmse, sqrt(res$mse))
})

# ---------------------------------------------------------------------------
# 7. sort_by_chr_and_start
# ---------------------------------------------------------------------------

test_that("sort_by_chr_and_start orders chromosomes correctly", {
  df <- data.frame(
    CHR   = c("chr2", "chr1", "chr10", "chrX"),
    START = c(100, 200, 50, 300),
    VAL   = 1:4
  )
  sorted <- SEMseeker:::sort_by_chr_and_start(df)
  expect_equal(as.character(sorted$CHR[1]), "chr1")
  expect_equal(as.character(sorted$CHR[2]), "chr2")
  expect_equal(as.character(sorted$CHR[3]), "chr10")
})

test_that("sort_by_chr_and_start sorts by START within same chromosome", {
  df <- data.frame(
    CHR   = c("chr1", "chr1", "chr1"),
    START = c(300, 100, 200),
    VAL   = 1:3
  )
  sorted <- SEMseeker:::sort_by_chr_and_start(df)
  expect_equal(sorted$START, c(100, 200, 300))
})

test_that("sort_by_chr_and_start handles empty data.frame", {
  df <- data.frame(CHR = character(0), START = integer(0))
  result <- SEMseeker:::sort_by_chr_and_start(df)
  expect_equal(nrow(result), 0)
})

# ---------------------------------------------------------------------------
# 8. compute_mean_delta_permutation_cpu
# ---------------------------------------------------------------------------

test_that("compute_mean_delta_permutation_cpu returns a numeric scalar", {
  df <- data.frame(
    burden = c(1.0, 1.2, 1.1, 2.0, 2.2, 2.1),
    group  = c(0, 0, 0, 1, 1, 1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::compute_mean_delta_permutation_cpu(f, df, shuffle = FALSE)
  expect_length(res, 1)
  expect_true(is.numeric(res))
})

test_that("compute_mean_delta_permutation_cpu: group 1 > group 0 gives negative delta", {
  # group 1 has higher values → mean(g1) - mean(g0) > 0 → function returns -1 * diff
  df <- data.frame(
    burden = c(1, 1, 1, 3, 3, 3),
    group  = c(0, 0, 0, 1, 1, 1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::compute_mean_delta_permutation_cpu(f, df, shuffle = FALSE)
  expect_lt(res, 0)    # -1 * (mean(g1) - mean(g0)) = -1 * 2 = -2
})

test_that("compute_mean_delta_permutation_cpu: shuffle changes the result", {
  set.seed(123)
  df <- data.frame(
    burden = c(1, 1, 1, 10, 10, 10),
    group  = c(0, 0, 0,  1,  1,  1)
  )
  f      <- stats::as.formula("burden ~ group")
  obs    <- SEMseeker:::compute_mean_delta_permutation_cpu(f, df, shuffle = FALSE)
  perms  <- replicate(50, SEMseeker:::compute_mean_delta_permutation_cpu(f, df, shuffle = TRUE))
  # Shuffled results should not all equal the observed statistic
  expect_false(all(perms == obs))
})

# ---------------------------------------------------------------------------
# 9. compute_spearman_permutation
# ---------------------------------------------------------------------------

test_that("compute_spearman_permutation returns a numeric scalar in [-1, 1]", {
  df <- data.frame(
    burden = 1:10,
    group  = 1:10
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::compute_spearman_permutation(f, df, shuffle = FALSE)
  expect_length(res, 1)
  expect_true(is.numeric(res))
  expect_gte(res, -1)
  expect_lte(res,  1)
})

test_that("compute_spearman_permutation: perfectly correlated data gives rho = 1", {
  df <- data.frame(
    burden = 1:20,
    group  = 1:20
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::compute_spearman_permutation(f, df, shuffle = FALSE)
  expect_equal(unname(res), 1, tolerance = 1e-10)
})

test_that("compute_spearman_permutation: shuffle changes the result", {
  set.seed(42)
  df <- data.frame(
    burden = 1:30,
    group  = 1:30
  )
  f     <- stats::as.formula("burden ~ group")
  obs   <- SEMseeker:::compute_spearman_permutation(f, df, shuffle = FALSE)
  perms <- replicate(30, SEMseeker:::compute_spearman_permutation(f, df, shuffle = TRUE))
  expect_false(all(perms == obs))
})

# ---------------------------------------------------------------------------
# 10. compute_quantile_delta_permutation
# ---------------------------------------------------------------------------

test_that("compute_quantile_delta_permutation returns a numeric scalar", {
  df <- data.frame(
    burden = c(1, 2, 3, 10, 11, 12),
    group  = c(0, 0, 0,  1,  1,  1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.5)
  expect_length(res, 1)
  expect_true(is.numeric(res))
})

test_that("compute_quantile_delta_permutation: Q50 of group 1 > Q50 of group 0", {
  df <- data.frame(
    burden = c(1, 2, 3, 10, 11, 12),
    group  = c(0, 0, 0,  1,  1,  1)
  )
  f   <- stats::as.formula("burden ~ group")
  res <- SEMseeker:::compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.5)
  expect_gt(res, 0)   # median(g1) - median(g0) = 11 - 2 = 9 > 0
})

test_that("compute_quantile_delta_permutation: Q10 and Q90 give different results", {
  set.seed(1)
  df <- data.frame(
    burden = c(stats::rnorm(20, 0), stats::rnorm(20, 5)),
    group  = c(rep(0, 20), rep(1, 20))
  )
  f   <- stats::as.formula("burden ~ group")
  r10 <- SEMseeker:::compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.10)
  r90 <- SEMseeker:::compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.90)
  expect_false(isTRUE(all.equal(r10, r90)))
})

test_that("compute_quantile_delta_permutation: shuffle changes the result", {
  set.seed(99)
  df <- data.frame(
    burden = c(1, 2, 3, 100, 101, 102),
    group  = c(0, 0, 0,   1,   1,   1)
  )
  f     <- stats::as.formula("burden ~ group")
  obs   <- SEMseeker:::compute_quantile_delta_permutation(f, df, shuffle = FALSE, quantile = 0.5)
  perms <- replicate(30, SEMseeker:::compute_quantile_delta_permutation(f, df, shuffle = TRUE, quantile = 0.5))
  expect_false(all(perms == obs))
})

# ---------------------------------------------------------------------------
# 11. describe_dataframe
# ---------------------------------------------------------------------------

test_that("describe_dataframe: returns one row per column", {
  df <- data.frame(x = 1:5, y = letters[1:5], z = c(1.1, NA, 2.2, NA, 3.3))
  res <- SEMseeker:::describe_dataframe(df)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
})

test_that("describe_dataframe: Variable column matches input column names", {
  df <- data.frame(alpha = 1:3, beta = c("a", "b", "c"))
  res <- SEMseeker:::describe_dataframe(df)
  expect_equal(as.character(res$Variable), c("alpha", "beta"))
})

test_that("describe_dataframe: Missing_Values counts NAs correctly", {
  df <- data.frame(x = c(1, NA, 3, NA, 5))
  res <- SEMseeker:::describe_dataframe(df)
  expect_equal(res$Missing_Values, 2)
})

test_that("describe_dataframe: Missing_Values_Percent is 100*NA/n", {
  df <- data.frame(x = c(NA, NA, 1, 1))  # 50% missing
  res <- SEMseeker:::describe_dataframe(df)
  expect_equal(res$Missing_Values_Percent, 50)
})

test_that("describe_dataframe: Mean/Median/Min/Max are NA for non-numeric columns", {
  df <- data.frame(label = c("a", "b", "c"))
  res <- SEMseeker:::describe_dataframe(df)
  expect_true(is.na(res$Mean))
  expect_true(is.na(res$Median))
  expect_true(is.na(res$Min))
  expect_true(is.na(res$Max))
})

test_that("describe_dataframe: numeric statistics are correct", {
  df <- data.frame(v = c(1, 2, 3, 4, 5))
  res <- SEMseeker:::describe_dataframe(df)
  expect_equal(res$Mean,   3)
  expect_equal(res$Median, 3)
  expect_equal(res$Min,    1)
  expect_equal(res$Max,    5)
})

# ---------------------------------------------------------------------------
# 12. polynomial_formula_build
# ---------------------------------------------------------------------------

test_that("polynomial_formula_build: degree-1 without covariates gives y ~ I(x^1)", {
  f <- SEMseeker:::polynomial_formula_build("y", "x", 1, character(0))
  expect_equal(deparse(f), "y ~ I(x^1)")
})

test_that("polynomial_formula_build: degree-3 without covariates has 3 terms", {
  f    <- SEMseeker:::polynomial_formula_build("y", "x", 3, character(0))
  fstr <- paste(deparse(f), collapse = " ")
  expect_true(grepl("I(x^1)", fstr, fixed = TRUE))
  expect_true(grepl("I(x^2)", fstr, fixed = TRUE))
  expect_true(grepl("I(x^3)", fstr, fixed = TRUE))
})

test_that("polynomial_formula_build: with covariates includes interaction terms", {
  f    <- SEMseeker:::polynomial_formula_build("y", "x", 2, c("age", "sex"))
  fstr <- paste(deparse(f), collapse = " ")  # deparse may split long formulas
  expect_true(grepl("I(x^1):age", fstr, fixed = TRUE))
  expect_true(grepl("I(x^2):sex", fstr, fixed = TRUE))
})

test_that("polynomial_formula_build: returns a formula object", {
  f <- SEMseeker:::polynomial_formula_build("outcome", "predictor", 2, character(0))
  expect_s3_class(f, "formula")
})

test_that("polynomial_formula_build: LHS is the dependent variable", {
  f <- SEMseeker:::polynomial_formula_build("response", "x", 2, character(0))
  expect_equal(as.character(f[[2]]), "response")
})
