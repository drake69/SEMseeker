test_that("coverage-analysis", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  # coverage_analysis strips POSITION and PROBE from the keys, so init_env must
  # receive at least one other area (GENE, ISLAND, …) or it returns early with 0 rows.
  semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"))

  tt <- semseeker:::get_meth_tech(signal_data)
  res_cov <- semseeker:::coverage_analysis(signal_data)
  testthat::expect_true(nrow(res_cov)!=0)

  # # test cases of some are fully missed
  # nitem <- 1e2
  # nsamples <- 21
  #
  # probe_features <- semseeker::pp_tot
  # probe_features <- probe_features[!is.na(probe_features$START),c("CHR","START","PROBE")]
  # probe_features <- unique(probe_features)
  # probe_features$END <- probe_features$START
  #
  # nitem <- min(nitem, nrow(probe_features))
  # probe_features <- probe_features[probe_features$PROBE %in% sample(x=probe_features[,"PROBE"] , size=nitem),]
  #
  # signal_data <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
  # signal_data <- as.data.frame(matrix(signal_data,nitem,nsamples))
  #
  # row.names(signal_data) <- probe_features$PROBE
  # ####################################################################################
  #
  # res_cov <- semseeker:::coverage_analysis(signal_data)
  # testthat::expect_true(nrow(res_cov)!=0)

  ####################################################################################

  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})

