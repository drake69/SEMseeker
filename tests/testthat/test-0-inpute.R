test_that("inpute", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]

  # Build a local copy with controlled NAs so this test is independent of whether
  # real or synthetic data is used.  We inject:
  #   - partial NAs in 3 rows (1-2 NAs each, well below the 10% removal threshold)
  #     → these are imputed in-place, row count unchanged
  #   - all-NA in 2 rows → removed (also caught by the >10% threshold removal)
  # With nsamples=30 the threshold is 0.1*30=3, so rows with ≤2 NAs are safe.
  local_sig <- signal_data[1:20, ]
  local_sig[1,  1]    <- NA   # 1 NA  → imputed
  local_sig[5,  1:2]  <- NA   # 2 NAs → imputed
  local_sig[10, 1]    <- NA   # 1 NA  → imputed
  local_sig[15, ]     <- NA   # all-NA row → removed
  local_sig[18, ]     <- NA   # all-NA row → removed
  nrow_local   <- nrow(local_sig)
  nrow_missed  <- sum(apply(local_sig, 1, function(x) all(is.na(x))))  # = 2
  ex_ante      <- sum(is.na(local_sig))

  ####################################################################################

  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute="median")
  signal_data_median <- SEMseeker:::inpute_missing_values(local_sig)

  nrows_ex_post <- nrow(signal_data_median)
  ex_post <- sum(is.na(signal_data_median))

  testthat::expect_true(ex_ante != 0)
  testthat::expect_true(ex_post == 0)
  testthat::expect_true(nrow_local - nrows_ex_post == nrow_missed)

  ####################################################################################

  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, iqrTimes = iqrTimes, inpute="mean")
  signal_data_mean <- SEMseeker:::inpute_missing_values(local_sig)

  nrows_ex_post <- nrow(signal_data_mean)
  ex_post <- sum(is.na(signal_data_mean))

  testthat::expect_true(ex_ante != 0)
  testthat::expect_true(ex_post == 0)
  testthat::expect_true(nrow_local - nrows_ex_post == nrow_missed)

  ####################################################################################

  SEMseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

