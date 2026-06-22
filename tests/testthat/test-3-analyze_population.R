test_that("analyze_population", {

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  ####################################################################################

  tt <- SEMseeker:::get_meth_tech(signal_data)

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- SEMseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  mySampleSheet <- mySampleSheet[mySampleSheet$Sample_Group!="Reference",]
  sp <- SEMseeker:::analyze_population(
    signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet,
    probe_features = probe_features
  )

  # analyze_population writes bed files rather than returning a data.frame; verify output was created
  data_dir <- file.path(tempFolder, "Data")
  testthat::expect_true(dir.exists(data_dir))
  testthat::expect_true(length(list.files(data_dir, recursive = TRUE)) > 0)


  ####################################################################################

  SEMseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)

})

