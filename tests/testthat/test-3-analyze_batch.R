test_that("analize_batch", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder,
    parallel_strategy = parallel_strategy,
    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    inpute="median"
  )

  ####################################################################################
  tt <- semseeker:::get_meth_tech(signal_data)
  ####################################################################################

  batch_id <- 2
  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  # ss <- mySampleSheet[mySampleSheet$Sample_Group!="Reference",]

  sp <- semseeker:::analyze_batch(
    signal_data =  signal_data,
    sample_sheet =  mySampleSheet
  )

  # analyze_batch writes files rather than returning a data.frame; verify output was created
  data_dir <- file.path(tempFolder, "Data")
  testthat::expect_true(dir.exists(data_dir))
  testthat::expect_true(length(list.files(data_dir, recursive = TRUE)) > 0)

  ####################################################################################
  semseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

