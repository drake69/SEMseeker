test_that("deltar_single_sample",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(result_folder= tempFolder, bonferroni_threshold = 0.05, showprogress=showprogress, start_fresh = TRUE, inpute="median")

  ####################################################################################

  gg <- SEMseeker:::get_meth_tech(signal_data)
  if (!exists("signal_thresholds"))
  {
    signal_data <- SEMseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- SEMseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]
  ####################################################################################
  # for (s in 1:ncol(signal_data))
  {
    s <- 2
    message(s)
    sample_id <- colnames(signal_data)[s]
    # Look up the actual sample group (don't hardcode to "Reference" — sample 2 may not be Reference)
    sample_detail <- mySampleSheet[mySampleSheet$Sample_ID == sample_id, c("Sample_ID","Sample_Group")]

    # Build a bed-like data.frame matching the production code path (values read from bed file)
    values_df <- data.frame(
      CHR   = probe_features$CHR[match(rownames(signal_data), probe_features$PROBE)],
      START = probe_features$START[match(rownames(signal_data), probe_features$PROBE)],
      END   = probe_features$END[match(rownames(signal_data), probe_features$PROBE)],
      VALUE = as.numeric(signal_data[, sample_id])
    )

    dr1 <- SEMseeker:::deltar_single_sample(
      values = values_df,
      thresholds = signal_thresholds,
      sample_detail = sample_detail
    )

    ds1 <- SEMseeker:::delta_single_sample(
      values = values_df,
      thresholds = signal_thresholds,
      sample_detail = sample_detail
    )

    ass <- SEMseeker:::analyze_single_sample(
      values = values_df,
      thresholds = signal_thresholds,
      figure = "HYPER",
      sample_detail = sample_detail
    )

    SEMseeker:::analyze_single_sample(
      values = values_df,
      thresholds = signal_thresholds,
      figure = "HYPO",
      sample_detail = sample_detail
    )

    mb <- SEMseeker:::analyze_single_sample_both(sample_detail,"MUTATIONS")

    result_folderData  <-  SEMseeker:::dir_check_and_create(tempFolder, "Data")
    outputFolder <- SEMseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"MUTATIONS_HYPER"))
    fileName <- SEMseeker:::file_path_build(outputFolder,c(sample_id,"MUTATIONS","HYPER"), "bed", add_gz = TRUE)
    mutations <- read.table(fileName)


    result_folderData  <-  SEMseeker:::dir_check_and_create(tempFolder, "Data")
    outputFolder <- SEMseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAR_HYPER"))
    fileName <- SEMseeker:::file_path_build(outputFolder,c(sample_id,"DELTAR","HYPER"), "bedgraph", add_gz = TRUE)
    deltar <- read.table(gzfile(fileName))
    # deltar <- read.csv(fileName, sep="\t")
    testthat::expect_true(nrow(deltar)>0)
    testthat::expect_true(nrow(deltar)==nrow(mutations))

    ####################################################################################

    outputFolder <- SEMseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAS_HYPER"))
    fileName <- SEMseeker:::file_path_build(outputFolder,c(sample_id,"DELTAS","HYPER"), "bedgraph", add_gz = TRUE)
    deltas <- read.table(fileName)
    testthat::expect_true(nrow(deltas)>0)
    testthat::expect_true(nrow(deltar)==nrow(deltas))

    ####################################################################################

    testthat::expect_true(nrow(mutations)==nrow(deltas))

    ####################################################################################

    outputFolder <- SEMseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAR_HYPER"))
    fileName <- SEMseeker:::file_path_build(outputFolder,c(sample_id,"DELTAR","HYPER"), "bedgraph", add_gz = TRUE)
    deltar <- read.table(fileName)
    testthat::expect_true(nrow(deltar)>0)

    outputFolder <- SEMseeker:::dir_check_and_create(result_folderData,c(sample_detail$Sample_Group,"DELTAS_HYPER"))
    fileName <- SEMseeker:::file_path_build(outputFolder,c(sample_id,"DELTAS","HYPER"), "bedgraph", add_gz = TRUE)
    deltas <- read.table(fileName)

    testthat::expect_true(nrow(deltar)==nrow(deltas))

  }

  ####################################################################################
  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)

})
