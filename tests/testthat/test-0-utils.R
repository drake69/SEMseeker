# Tests for utility helper functions (no session required)
#
# Covered:
#  - boolean_check()          string/logical/NA → logical coercion
#  - file_path_build()        path construction helper
#  - name_composer()          filename fragment assembler
#  - data.frame_add.column()  safe column add/update on data.frames
#  - filter_sql()             SQL-condition row filter (requires sqldf)

# ---------------------------------------------------------------------------
# 1. boolean_check
# ---------------------------------------------------------------------------

test_that("boolean_check returns FALSE for NA, empty string, NULL", {
  expect_false(SEMseeker:::boolean_check(NA))
  expect_false(SEMseeker:::boolean_check(""))
  expect_false(SEMseeker:::boolean_check(NULL))
})

test_that("boolean_check returns TRUE for TRUE, 'TRUE', 'T', 'true', 'True', '1', 1", {
  expect_true(SEMseeker:::boolean_check(TRUE))
  expect_true(SEMseeker:::boolean_check("TRUE"))
  expect_true(SEMseeker:::boolean_check("T"))
  expect_true(SEMseeker:::boolean_check("true"))
  expect_true(SEMseeker:::boolean_check("True"))
  expect_true(SEMseeker:::boolean_check("1"))
  expect_true(SEMseeker:::boolean_check(1))
})

test_that("boolean_check returns FALSE for FALSE, 'FALSE', '0', 0", {
  expect_false(SEMseeker:::boolean_check(FALSE))
  expect_false(SEMseeker:::boolean_check("FALSE"))
  expect_false(SEMseeker:::boolean_check("0"))
  expect_false(SEMseeker:::boolean_check(0))
})

test_that("boolean_check returns FALSE for arbitrary strings", {
  expect_false(SEMseeker:::boolean_check("yes"))
  expect_false(SEMseeker:::boolean_check("no"))
  expect_false(SEMseeker:::boolean_check("maybe"))
})

# ---------------------------------------------------------------------------
# 2. file_path_build
# ---------------------------------------------------------------------------

test_that("file_path_build assembles path correctly", {
  result <- SEMseeker:::file_path_build("/base", c("a", "b"), "csv")
  expect_true(grepl("^/base", result))
  expect_true(grepl("A_B\\.csv$", result))
})

test_that("file_path_build adds .gz when add_gz = TRUE", {
  result <- SEMseeker:::file_path_build("/base", c("sample"), "bed", add_gz = TRUE)
  expect_true(endsWith(result, ".gz"))
  expect_true(grepl("\\.bed\\.gz$", result))
})

test_that("file_path_build does not add .gz when add_gz = FALSE", {
  result <- SEMseeker:::file_path_build("/base", c("sample"), "csv", add_gz = FALSE)
  expect_false(endsWith(result, ".gz"))
})

test_that("file_path_build cleans the filename part (upper, underscore)", {
  result <- SEMseeker:::file_path_build("/tmp", c("hello world", "foo.bar"), "txt")
  fname <- basename(result)
  expect_true(grepl("HELLO_WORLD", fname))
  expect_true(grepl("FOO_BAR",     fname))
})

test_that("file_path_build collapses multiple filenames with underscore", {
  result <- SEMseeker:::file_path_build("/tmp", c("A", "B", "C"), "csv")
  expect_true(grepl("A_B_C", basename(result)))
})

test_that("file_path_build avoids double dots", {
  result <- SEMseeker:::file_path_build("/tmp", c("X"), "csv")
  expect_false(grepl("\\.\\.", result))
})

# ---------------------------------------------------------------------------
# 3. name_composer
# ---------------------------------------------------------------------------

test_that("name_composer joins fragments with underscore", {
  result <- SEMseeker:::name_composer("a", "b", "c")
  expect_equal(result, "a_b_c")
})

test_that("name_composer replaces spaces with underscore", {
  result <- SEMseeker:::name_composer("hello world", "foo")
  expect_false(grepl(" ", result))
  expect_true(grepl("_", result))
})

test_that("name_composer replaces colons with underscore", {
  result <- SEMseeker:::name_composer("key:value")
  expect_false(grepl(":", result))
})

test_that("name_composer collapses double underscores", {
  result <- SEMseeker:::name_composer("a", "", "b")
  expect_false(grepl("__", result))
})

test_that("name_composer avoids double dots", {
  result <- SEMseeker:::name_composer("file.name", ".csv")
  expect_false(grepl("\\.\\.", result))
})

test_that("name_composer handles a single fragment", {
  result <- SEMseeker:::name_composer("single")
  expect_equal(result, "single")
})

# ---------------------------------------------------------------------------
# 4. data.frame_add.column
# ---------------------------------------------------------------------------

test_that("data.frame_add.column adds a new column when absent", {
  df <- data.frame(x = 1:3)
  result <- SEMseeker:::data.frame_add.column(df, "y", c(4, 5, 6))
  expect_true("y" %in% colnames(result))
  expect_equal(result$y, c(4, 5, 6))
})

test_that("data.frame_add.column updates existing column", {
  df <- data.frame(x = 1:3, y = c(0, 0, 0))
  result <- SEMseeker:::data.frame_add.column(df, "y", c(7, 8, 9))
  expect_equal(result$y, c(7, 8, 9))
  expect_equal(ncol(result), 2)   # no extra column added
})

test_that("data.frame_add.column works on empty data.frame", {
  df <- data.frame()
  result <- SEMseeker:::data.frame_add.column(df, "col", 42)
  expect_true("col" %in% colnames(result))
})

test_that("data.frame_add.column preserves existing columns", {
  df <- data.frame(a = 1:3, b = 4:6)
  result <- SEMseeker:::data.frame_add.column(df, "c", 7:9)
  expect_true(all(c("a", "b", "c") %in% colnames(result)))
  expect_equal(result$a, 1:3)
})

# ---------------------------------------------------------------------------
# 5. filter_sql
# ---------------------------------------------------------------------------
# filter_sql uses sqldf internally and calls log_event (requires a session).
# Tests that exercise the sqldf path are wrapped with init_env / close_env.
# sqldf resolves table names from the calling frame; we therefore call
# filter_sql via a thin wrapper so that the local variable name matches
# what the function expects.

test_that("filter_sql returns data unchanged for empty conditions", {
  df <- data.frame(x = 1:5, g = c("a","b","a","b","a"))
  result <- SEMseeker:::filter_sql(c(), df)
  expect_equal(nrow(result), nrow(df))
})

test_that("filter_sql filters rows by numeric condition", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(x = 1:10)
  result <- SEMseeker:::filter_sql("x > 5", df)
  expect_true(all(result$x > 5))
  expect_equal(nrow(result), 5)

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("filter_sql filters rows by character condition", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(
    val = c(1, 2, 3),
    grp = c("case", "control", "case"),
    stringsAsFactors = FALSE
  )
  result <- SEMseeker:::filter_sql("grp = 'case'", df)
  expect_equal(nrow(result), 2)
  expect_true(all(result$grp == "case"))

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("filter_sql applies multiple conditions sequentially", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(x = 1:10, y = c(rep("A",5), rep("B",5)),
                   stringsAsFactors = FALSE)
  result <- SEMseeker:::filter_sql(c("x > 3", "x < 8"), df)
  expect_true(all(result$x > 3 & result$x < 8))

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("filter_sql returns empty data.frame when no rows match", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  df <- data.frame(x = 1:5)
  result <- SEMseeker:::filter_sql("x > 100", df)
  expect_equal(nrow(result), 0)

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
