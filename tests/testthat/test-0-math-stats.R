## test-0-math-stats.R
## Unit tests for pure mathematical / statistical utility functions:
##   ssim, variation_of_information, io_convertTextToNumeric,
##   util_split_and_clean, metrics_filter, normalize_minimize, normalize_maximize,
##   util_substitute_infinite, calculate_collinearity_score
##
## These functions have no file I/O or network dependencies and can be
## exercised directly after devtools::load_all().

# -----------------------------------------------------------------------
# ssim
# -----------------------------------------------------------------------
test_that("ssim: identical vectors return value ~1", {
  x <- c(1, 2, 3, 4, 5)
  result <- SEMseeker:::ssim(x, x)
  expect_true(result > 0.999)
})

test_that("ssim: identical matrices return value ~1", {
  m <- matrix(1:9, nrow = 3)
  result <- SEMseeker:::ssim(m, m)
  expect_true(result > 0.999)
})

test_that("ssim: result is in a reasonable range for similar vectors", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1.1, 2.1, 2.9, 4.0, 5.1)
  result <- SEMseeker:::ssim(x, y)
  expect_true(result > 0.9)
  expect_true(result <= 1.0)
})

test_that("ssim: very different vectors return lower value than similar ones", {
  x <- c(1, 2, 3, 4, 5)
  y_similar  <- x + 0.1
  y_different <- rev(x) * 10
  ssim_similar  <- SEMseeker:::ssim(x, y_similar)
  ssim_different <- SEMseeker:::ssim(x, y_different)
  expect_true(ssim_similar > ssim_different)
})

test_that("ssim: mismatched lengths raise an error", {
  expect_error(SEMseeker:::ssim(1:3, 1:5))
})

test_that("ssim: mismatched matrix dimensions raise an error", {
  expect_error(SEMseeker:::ssim(matrix(1:4, 2), matrix(1:9, 3)))
})

# -----------------------------------------------------------------------
# variation_of_information
# -----------------------------------------------------------------------
test_that("variation_of_information: identical partitions return 0", {
  x <- c(1, 1, 2, 2, 3, 3)
  result <- SEMseeker:::variation_of_information(x, x)
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("variation_of_information: is non-negative", {
  x <- c(1, 1, 2, 2, 3, 3)
  y <- c(1, 2, 1, 3, 2, 3)
  result <- SEMseeker:::variation_of_information(x, y)
  expect_true(result >= 0)
})

test_that("variation_of_information: different partitions return >0", {
  x <- c(1, 1, 2, 2)
  y <- c("a", "b", "a", "b")
  result <- SEMseeker:::variation_of_information(x, y)
  expect_true(result > 0)
})

test_that("variation_of_information: symmetric (VI(x,y) == VI(y,x))", {
  x <- c(1, 1, 2, 2, 3, 3)
  y <- c(1, 2, 2, 3, 3, 1)
  expect_equal(
    SEMseeker:::variation_of_information(x, y),
    SEMseeker:::variation_of_information(y, x),
    tolerance = 1e-10
  )
})

test_that("variation_of_information: perfectly correlated partitions give 0", {
  x <- c("a", "a", "b", "b")
  y <- c(1,   1,   2,   2  )
  result <- SEMseeker:::variation_of_information(x, y)
  expect_equal(result, 0, tolerance = 1e-10)
})

# -----------------------------------------------------------------------
# io_convertTextToNumeric
# -----------------------------------------------------------------------
test_that("io_convertTextToNumeric: plain decimal with dot", {
  expect_equal(SEMseeker:::io_convertTextToNumeric("3.14"), 3.14)
})

test_that("io_convertTextToNumeric: European decimal comma", {
  expect_equal(SEMseeker:::io_convertTextToNumeric("3,14"), 3.14)
})

test_that("io_convertTextToNumeric: integer string", {
  expect_equal(SEMseeker:::io_convertTextToNumeric("42"), 42)
})

test_that("io_convertTextToNumeric: US thousands separator (1,200,000.34)", {
  result <- SEMseeker:::io_convertTextToNumeric("1,200,000.34")
  expect_equal(result, 1200000.34)
})

test_that("io_convertTextToNumeric: European thousands separator (1.200.000,34)", {
  result <- SEMseeker:::io_convertTextToNumeric("1.200.000,34")
  expect_equal(result, 1200000.34)
})

# -----------------------------------------------------------------------
# util_split_and_clean
# -----------------------------------------------------------------------
test_that("util_split_and_clean: splits on default '+' delimiter", {
  result <- SEMseeker:::util_split_and_clean("a+b+c")
  expect_equal(sort(result), c("a", "b", "c"))
})

test_that("util_split_and_clean: single element returns itself", {
  expect_equal(SEMseeker:::util_split_and_clean("hello"), "hello")
})

test_that("util_split_and_clean: trims whitespace around parts", {
  result <- SEMseeker:::util_split_and_clean("a + b + c")
  expect_equal(sort(result), c("a", "b", "c"))
})

test_that("util_split_and_clean: removes duplicates", {
  result <- SEMseeker:::util_split_and_clean("x+x+y")
  expect_equal(sort(result), c("x", "y"))
})

test_that("util_split_and_clean: custom split character", {
  result <- SEMseeker:::util_split_and_clean("a,b,c", split = ",")
  expect_equal(sort(result), c("a", "b", "c"))
})

test_that("util_split_and_clean: empty string returns character(0) after whitespace filter", {
  result <- SEMseeker:::util_split_and_clean("")
  # "" is not NULL/NA-empty (length 1), so rlang::is_empty is FALSE;
  # trimws("") == "" is then removed → character(0)
  expect_equal(length(result), 0)
})

# -----------------------------------------------------------------------
# metrics_filter
# -----------------------------------------------------------------------
test_that("metrics_filter: 'none' transformation returns all metrics unchanged", {
  all_metrics <- c("MAE", "RMSE", "COUNT_SIGN", "R_SQUARED")
  result <- SEMseeker:::metrics_filter(all_metrics, "none")
  expect_equal(sort(result), sort(all_metrics))
})

test_that("metrics_filter: 'scale' transformation returns all metrics (sorted unique)", {
  all_metrics <- c("MAE", "COUNT_SIGN", "R_SQUARED", "MAE")  # duplicate
  result <- SEMseeker:::metrics_filter(all_metrics, "scale")
  expect_true(length(result) <= length(all_metrics))
  expect_true("MAE" %in% result)
})

test_that("metrics_filter: 'log' removes scale-affected metrics", {
  affected <- toupper(SEMseeker::metrics_properties[
    SEMseeker::metrics_properties$Affected_by_Scaling == TRUE, "Metric"])
  not_affected <- toupper(SEMseeker::metrics_properties[
    SEMseeker::metrics_properties$Affected_by_Scaling == FALSE, "Metric"])

  all_metrics <- c(affected[1], not_affected[1])
  result <- SEMseeker:::metrics_filter(all_metrics, "log")

  expect_false(affected[1] %in% result)
  expect_true(not_affected[1] %in% result)
})

test_that("metrics_filter: 'none' returns exactly the original vector (no reorder)", {
  metrics <- c("RMSE", "MAE", "COUNT_SIGN")
  result <- SEMseeker:::metrics_filter(metrics, "none")
  expect_equal(result, metrics)
})

test_that("metrics_filter: non-'none' transformations return a sorted unique result", {
  metrics <- c("RMSE", "MAE", "COUNT_SIGN", "MAE")   # duplicate
  result <- SEMseeker:::metrics_filter(metrics, "log")
  expect_equal(result, sort(unique(result)))
})

# -----------------------------------------------------------------------
# normalize_minimize / normalize_maximize
# -----------------------------------------------------------------------
test_that("normalize_minimize: min value maps to 1, max to 0", {
  x <- c(1, 2, 3, 4, 5)
  r <- SEMseeker:::normalize_minimize(x)
  expect_equal(r[which.min(x)], 1)
  expect_equal(r[which.max(x)], 0)
})

test_that("normalize_minimize: all values in [0, 1]", {
  x <- c(10, 20, 15, 5, 25)
  r <- SEMseeker:::normalize_minimize(x)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("normalize_minimize: range is [0, 1]", {
  x <- c(3, 7, 1, 9, 5)
  r <- SEMseeker:::normalize_minimize(x)
  expect_equal(min(r), 0)
  expect_equal(max(r), 1)
})

test_that("normalize_minimize: monotonically decreasing with input", {
  x <- 1:10
  r <- SEMseeker:::normalize_minimize(x)
  expect_true(all(diff(r) < 0))
})

test_that("normalize_maximize: min value maps to 0, max to 1", {
  x <- c(1, 2, 3, 4, 5)
  r <- SEMseeker:::normalize_maximize(x)
  expect_equal(r[which.min(x)], 0)
  expect_equal(r[which.max(x)], 1)
})

test_that("normalize_maximize: all values in [0, 1]", {
  x <- c(10, 20, 15, 5, 25)
  r <- SEMseeker:::normalize_maximize(x)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("normalize_maximize: range is [0, 1]", {
  x <- c(3, 7, 1, 9, 5)
  r <- SEMseeker:::normalize_maximize(x)
  expect_equal(min(r), 0)
  expect_equal(max(r), 1)
})

test_that("normalize_maximize: monotonically increasing with input", {
  x <- 1:10
  r <- SEMseeker:::normalize_maximize(x)
  expect_true(all(diff(r) > 0))
})

test_that("normalize_minimize and normalize_maximize are complementary", {
  x <- c(2, 5, 1, 8, 3)
  r_min <- SEMseeker:::normalize_minimize(x)
  r_max <- SEMseeker:::normalize_maximize(x)
  expect_equal(r_min + r_max, rep(1, length(x)))
})

# -----------------------------------------------------------------------
# util_substitute_infinite  (calls log_event — uses global tempFolders from setup.R)
# -----------------------------------------------------------------------
test_that("util_substitute_infinite: no Inf values returns unchanged data frame", {
  SEMseeker:::init_env(result_folder = tempFolders[1], start_fresh = TRUE)
  df <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
  result <- SEMseeker:::util_substitute_infinite(df)
  expect_equal(as.numeric(result$x), c(1, 2, 3))
  expect_equal(as.numeric(result$y), c(4, 5, 6))
  SEMseeker:::close_env()
})

test_that("util_substitute_infinite: +Inf replaced by max finite value", {
  SEMseeker:::init_env(result_folder = tempFolders[2], start_fresh = TRUE)
  df <- data.frame(x = c(1, Inf, 3))
  result <- SEMseeker:::util_substitute_infinite(df)
  expect_false(any(is.infinite(as.numeric(result$x))))
  expect_equal(as.numeric(result$x[2]), 3)  # max finite = 3, sign(+Inf)=+1
  SEMseeker:::close_env()
})

test_that("util_substitute_infinite: -Inf replaced by negative max finite value", {
  SEMseeker:::init_env(result_folder = tempFolders[3], start_fresh = TRUE)
  df <- data.frame(x = c(1, -Inf, 3))
  result <- SEMseeker:::util_substitute_infinite(df)
  expect_false(any(is.infinite(as.numeric(result$x))))
  expect_equal(as.numeric(result$x[2]), -3)  # max finite = 3, sign(-Inf)=-1
  SEMseeker:::close_env()
})

test_that("util_substitute_infinite: returns a data.frame", {
  SEMseeker:::init_env(result_folder = tempFolders[4], start_fresh = TRUE)
  df <- data.frame(x = c(1, Inf, 2), y = c(-Inf, 3, 4))
  result <- SEMseeker:::util_substitute_infinite(df)
  expect_s3_class(result, "data.frame")
  SEMseeker:::close_env()
})

# -----------------------------------------------------------------------
# calculate_collinearity_score  (calls log_event — uses tempFolders from setup.R)
# -----------------------------------------------------------------------
test_that("calculate_collinearity_score: uncorrelated variables → no removal", {
  set.seed(42)
  SEMseeker:::init_env(result_folder = tempFolders[5], start_fresh = TRUE)
  df <- data.frame(
    a = stats::rnorm(50),
    b = stats::rnorm(50),
    c = stats::rnorm(50)
  )
  result <- SEMseeker:::calculate_collinearity_score(df)
  expect_equal(length(result), 0)
  SEMseeker:::close_env()
})

test_that("calculate_collinearity_score: highly correlated variables flagged for removal", {
  set.seed(42)
  SEMseeker:::init_env(result_folder = tempFolders[6], start_fresh = TRUE)
  base <- stats::rnorm(50)
  df <- data.frame(
    a = base,
    b = base + stats::rnorm(50, sd = 0.01),   # almost identical to a
    c = stats::rnorm(50)
  )
  result <- SEMseeker:::calculate_collinearity_score(df)
  expect_true(length(result) > 0)
  SEMseeker:::close_env()
})

test_that("calculate_collinearity_score: returns empty vector when no collinearity", {
  set.seed(1)
  SEMseeker:::init_env(result_folder = tempFolders[7], start_fresh = TRUE)
  df <- data.frame(x = stats::rnorm(30), y = stats::rnorm(30))
  result <- SEMseeker:::calculate_collinearity_score(df)
  # No collinearity → returns c() (zero-length vector)
  expect_equal(length(result), 0)
  SEMseeker:::close_env()
})
