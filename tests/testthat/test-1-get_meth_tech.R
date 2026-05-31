test_that("get-meth_tech", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  SEMseeker:::init_env(result_folder =  tempFolder, parallel_strategy = parallel_strategy, maxResources = 90, figures = "HYPER", markers = "DELTAS", areas = "GENE")

  ####################################################################################

  # signal_data_27 <- subset(signal_data, rownames(signal_data) %in% probe_features[probe_features$K27,"PROBE"])
  ssEnv <- SEMseeker:::get_meth_tech(signal_data)
  testthat::expect_true(ssEnv$tech!="K27")

  ####################################################################################

  # signal_data_450 <- subset(signal_data, rownames(signal_data) %in% probe_features[probe_features$K450,"PROBE"])
  ssEnv <- SEMseeker:::get_meth_tech(signal_data)
  # synthetic probes sampled from PROBES are predominantly K850 (EPIC array)
  testthat::expect_true(ssEnv$tech=="K850")

  ####################################################################################

  # signal_data_850 <- subset(signal_data, rownames(signal_data) %in% probe_features[probe_features$K850,"PROBE"])
  ssEnv <- SEMseeker:::get_meth_tech(signal_data)
  testthat::expect_true(ssEnv$tech!="K27")

  ####################################################################################
  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)

})
