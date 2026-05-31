# Tests for the public semseeker() dispatcher and its helpers.
# The dispatcher is covered by unit tests; the full end-to-end pipeline is
# exercised by test-6-semseeker.R.

# ---- .detect_input_type -------------------------------------------------

test_that(".detect_input_type recognises bedmethyl file paths", {
  expect_equal(semseeker:::.detect_input_type("sample.bed"),       "bedmethyl")
  expect_equal(semseeker:::.detect_input_type("sample.tsv"),       "bedmethyl")
  expect_equal(semseeker:::.detect_input_type("sample.bedmethyl"), "bedmethyl")
  expect_equal(semseeker:::.detect_input_type(c("a.bed","b.bed")), "bedmethyl")
})

test_that(".detect_input_type recognises coord data frames", {
  df <- data.frame(CHR = "chr1", START = 100L, s1 = 0.5, s2 = 0.7)
  expect_equal(semseeker:::.detect_input_type(df), "coord_df")
})

test_that(".detect_input_type recognises probe-indexed matrices", {
  m <- matrix(runif(6), nrow = 3, dimnames = list(
    c("cg00000029","cg00000165","cg00000236"), c("s1","s2")))
  expect_equal(semseeker:::.detect_input_type(m), "matrix")

  df <- data.frame(s1 = c(0.1, 0.2), s2 = c(0.3, 0.4),
                   row.names = c("cg00000029", "cg00000165"))
  expect_equal(semseeker:::.detect_input_type(df), "matrix")
})

test_that(".detect_input_type rejects unsupported extensions and classes", {
  expect_error(semseeker:::.detect_input_type("data.csv"),  "unsupported extensions")
  expect_error(semseeker:::.detect_input_type(list(1, 2)),  "cannot infer input_type")
})

# ---- .looks_like_mvalues ------------------------------------------------

test_that(".looks_like_mvalues is FALSE for beta values in [0,1]", {
  m <- matrix(runif(20), nrow = 5, dimnames = list(
    paste0("cg", sprintf("%08d", seq_len(5))), paste0("s", 1:4)))
  expect_false(semseeker:::.looks_like_mvalues(m))
})

test_that(".looks_like_mvalues is TRUE when abs values exceed 1 (M-values)", {
  m <- matrix(rnorm(20, mean = 0, sd = 3), nrow = 5, dimnames = list(
    paste0("cg", sprintf("%08d", seq_len(5))), paste0("s", 1:4)))
  expect_true(semseeker:::.looks_like_mvalues(m))
})

test_that(".looks_like_mvalues ignores CHR/START/END columns", {
  df <- data.frame(
    CHR = "chr1", START = 100L, END = 101L,
    s1 = c(0.1, 0.2, 0.3), s2 = c(0.4, 0.5, 0.6))
  expect_false(semseeker:::.looks_like_mvalues(df))
})

# ---- mvalue_to_beta -----------------------------------------------------

test_that("mvalue_to_beta converts a matrix via 2^M / (1 + 2^M)", {
  m <- matrix(c(-2, 0, 2), ncol = 1)
  beta <- semseeker:::mvalue_to_beta(m)
  expect_equal(as.numeric(beta), c(0.2, 0.5, 0.8), tolerance = 1e-6)
})

test_that("mvalue_to_beta preserves CHR/START/END columns in coord data frames", {
  df <- data.frame(CHR = "chr1", START = 100L, END = 101L,
                   s1 = c(-2, 0, 2),
                   stringsAsFactors = FALSE)
  out <- semseeker:::mvalue_to_beta(df, coord_cols = c("CHR","START","END"))
  expect_equal(out$CHR,   rep("chr1", 3))
  expect_equal(out$START, rep(100L, 3))
  expect_equal(out$s1,    c(0.2, 0.5, 0.8), tolerance = 1e-6)
})

# ---- .validate_tech_build -----------------------------------------------

test_that(".validate_tech_build stops on LONGREAD+hg19 when strict=TRUE", {
  expect_error(
    semseeker:::.validate_tech_build("LONGREAD", "hg19", strict = TRUE),
    "LONGREAD.*hg19"
  )
})

test_that(".validate_tech_build warns (not stops) when strict=FALSE", {
  expect_warning(
    semseeker:::.validate_tech_build("LONGREAD", "hg19", strict = FALSE),
    "LONGREAD.*hg19"
  )
})

test_that(".validate_tech_build accepts LONGREAD+hg38", {
  expect_silent(semseeker:::.validate_tech_build("LONGREAD", "hg38", strict = TRUE))
})

test_that(".validate_tech_build is a no-op when tech is NULL or empty", {
  expect_silent(semseeker:::.validate_tech_build(NULL, "hg19", strict = TRUE))
  expect_silent(semseeker:::.validate_tech_build("",   "hg19", strict = TRUE))
})

# ---- semseeker() surface: strict_build_check error path -----------------

test_that("semseeker() fails fast on LONGREAD+hg19 when strict_build_check=TRUE", {
  ss <- data.frame(Sample_ID = "s1", Sample_Group = "Case",
                   stringsAsFactors = FALSE)
  m  <- matrix(runif(4), nrow = 2,
               dimnames = list(c("cg00000029","cg00000165"), c("s1","s2")))
  expect_error(
    semseeker::semseeker(
      input = m,
      sample_sheet = ss,
      result_folder = tempfile(),
      tech = "LONGREAD",
      genome_build = "hg19",
      strict_build_check = TRUE
    ),
    "LONGREAD.*hg19"
  )
})
