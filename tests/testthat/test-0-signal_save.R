test_that("signal-save",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  # In the normal pipeline analyze_batch() calls get_meth_tech() before signal_save()
  # so that ssEnv$tech is already set when probe_features_get() runs inside signal_save().
  SEMseeker:::get_meth_tech(signal_data)
  SEMseeker:::signal_save(signal_data,mySampleSheet,batch_id )

  # signal_save writes the probe-level parquet with subarea "WHOLE"
  signal_file <- SEMseeker:::pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")
  testthat::expect_true(file.exists(signal_file))
  # it also writes the position-level file
  position_file <- SEMseeker:::pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION", "WHOLE")
  testthat::expect_true(file.exists(position_file))


})
