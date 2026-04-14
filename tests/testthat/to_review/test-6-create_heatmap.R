test_that("create_heatmap", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  ####################################################################################

  SEMseeker:::get_meth_tech(signal_data)

  ####################################################################################

  if (!exists("signal_thresholds"))
  {
    signal_data <- SEMseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  batch_id <- 1
  sp <- SEMseeker:::analyze_batch(signal_data =  signal_data,
    sample_sheet =  mySampleSheet,

    bonferroni_threshold =  bonferroni_threshold,
    iqrTimes =  iqrTimes,
    batch_id = batch_id
  )

  SEMseeker:::create_multiple_bed( sample_sheet = mySampleSheet)
  sp$Sample_Group <- mySampleSheet$Sample_Group


  SEMseeker:::create_excel_pivot()
  SEMseeker:::create_heatmap()

  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/CHR/Control_Vs_Case_CHR_MUTATIONS_HYPER.png")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_Vs_Case_GENE_AREA_MUTATIONS_HYPER.png")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE_AREA/Control_Vs_Case_GENE_AREA_DELTAS_HYPER.png")))
  testthat::expect_true(file.exists(file.path(ssEnv$result_folderChart,"/GENE/Control_Vs_Case_GENE_DELTAS_HYPER.png")))

  # unlink(tempFolder,recursive = TRUE)
  SEMseeker:::close_env()
})
