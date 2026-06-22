test_that("get-meth_tech", {

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  # start_fresh = TRUE ensures no stale session_info.rds from a previous test
  # (all 0-* tests share tempFolders[1] in devtools::test()) affects tech detection.
  SEMseeker:::init_env(result_folder = tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "HYPER", markers = "DELTAS", areas = "GENE", start_fresh = TRUE)

  ####################################################################################

  # signal_data_27 <- subset(signal_data, rownames(signal_data) %in% probe_features[probe_features$K27,"PROBE"])
  ssEnv <- SEMseeker:::get_meth_tech(signal_data)
  testthat::expect_true(ssEnv$tech!="K27")

  ####################################################################################

  # signal_data_450 <- subset(signal_data, rownames(signal_data) %in% probe_features[probe_features$K450,"PROBE"])
  ssEnv <- SEMseeker:::get_meth_tech(signal_data)
  # GSE133774 uses EPIC 850k probes → detect_tech_from_anno returns K850
  testthat::expect_true(ssEnv$tech=="K850")

  ####################################################################################

  # signal_data_850 <- subset(signal_data, rownames(signal_data) %in% probe_features[probe_features$K850,"PROBE"])
  ssEnv <- SEMseeker:::get_meth_tech(signal_data)
  testthat::expect_true(ssEnv$tech!="K27")

  ####################################################################################
  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)

})
