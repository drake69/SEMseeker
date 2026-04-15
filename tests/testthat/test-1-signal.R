test_that("signal",{

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, inpute="median")

  sample_detail <- mySampleSheet[1,]

  res <- SEMseeker:::signal_single_sample(values = signal_data[,sample_detail$Sample_ID], sample_detail = sample_detail, probe_features = probe_features)

  folder_to_save <- SEMseeker:::dir_check_and_create(ssEnv$result_folderData,c(sample_detail$Sample_Group,paste0("SIGNAL","_", "MEAN", sep = "")))
  signal_file <- SEMseeker:::file_path_build(folder_to_save,c(sample_detail$Sample_ID,"SIGNAL","MEAN"),"bedgraph", add_gz = TRUE)
  testthat::expect_true(file.exists(signal_file))

  signal_file <- read.table(signal_file, sep="\t")
  testthat::expect_true(nrow(signal_file)==nrow(signal_data))


})
