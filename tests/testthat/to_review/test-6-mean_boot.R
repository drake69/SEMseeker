test_that("association_analysis", {


  # tempFolder <- tempFolders[1]
  # tempFolders <- tempFolders[-1]
  #
  # ####################################################################################
  #
  # ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)
  # semseeker:::semseeker( sample_sheet =  mySampleSheet,signal_data =  signal_data,result_folder = tempFolder,parallel_strategy=parallel_strategy, figures="HYPER",
  #   markers=c("DELTAQ"), areas=c("PROBE","CHR","GENE"))
  #
  # ####################################################################################
  #
  # inferenceFolder <- file.path(tempFolder,"Inference")
  # inference_details <- expand.grid("independent_variable"= "Sample_Group",
  #   "family_test"=c("mean-bootstrap_100_1000_0.95"),
  #   "transformation"="",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="HYPER",
  #   markers=c("DELTAS","DELTAQ"), areas="GENE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)

  # ####################################################################################
  #
  # inferenceFolder <- file.path(tempFolder,"Inference")
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("poisson"),
  #   "transformation"="",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="HYPER",
  #   markers=c("DELTAS","DELTAQ"), areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_poisson_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inferenceFolder <- file.path(tempFolder,"Inference")
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5_5_10_0.9"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="HYPER",
  #   markers=c("DELTAS","DELTAQ"), areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_5_10_0.9_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Covariates1",
  #                                  "covariates"=c("Phenotest+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5_5_10"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  # # inference_details,result_folder, maxResources, parallel_strategy
  # # test with area selection
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER",
  #   markers=c("DELTAS","DELTAQ"), areas="CHR",areas_selection=c("chr1","chr2"))
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Covariates1_scale_quantreg_0.5_5_10_Phenotest_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # firstRun <- nrow(localFileRes)
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER",
  #   markers=c("DELTAS","DELTAQ"), areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Covariates1_scale_quantreg_0.5_5_10_Phenotest_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  # secondRun <- nrow(localFileRes)
  # testthat::expect_true(secondRun > firstRun)
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  # # inference_details,result_folder, maxResources, parallel_strategy
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers=c("DELTAS","DELTAQ"),
  #   areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # # inference_details <- expand.grid("independent_variable"= "Phenotest",
  # #   "covariates"=c("Covariates1+Covariates2"),
  # #   "family_test"=c("quantreg_0.5"),
  # #   "transformation"="scale",
  # #   "depth_analysis"=3,
  # #   "filter_p_value" = FALSE)
  # #
  # # # inference_details,result_folder, maxResources, parallel_strategy
  # # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,
  # #   figures="HYPER", markers=c("DELTAS","DELTAQ"), areas="PROBE")
  # #
  # # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_Covariates1_Covariates2", extension = "csv")
  # # localFileRes <- read.table(fileToRead, sep=";")
  # # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5_5_10_0.99"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers=c("DELTAS","DELTAQ"), areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_5_10_0.99_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("quantreg_0.5_5_10_0.99_np"),
  #   "transformation"="scale",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers=c("DELTAS","DELTAQ"), areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_5_10_0.99_np_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # #todo: test incremental association analysis
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #                                  "covariates"=c("Covariates1+Covariates2"),
  #                                  "family_test"=c("quantreg_0.5_5_10_0.9"),
  #                                  "transformation"="scale",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers=c("DELTAS","DELTAQ"), areas="CHR")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_quantreg_0.5_5_10_0.9_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # test_both <- nrow(localFileRes)
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Group",
  #                                  "covariates"="",
  #                                  "family_test"=c("wilcoxon"),
  #                                  "transformation"="log",
  #                                  "depth_analysis"=1,
  #                                  "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "1_Group_log_wilcoxon_", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  #
  # inference_details <- expand.grid("independent_variable"= "Group",
  #   "covariates"="",
  #   "family_test"=c("t.test"),
  #   "transformation"="log",
  #   "depth_analysis"=1,
  #   "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "1_Group_log_wilcoxon_", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Group",
  #   "covariates"="",
  #   "family_test"=c("pearson"),
  #   "transformation"="log",
  #   "depth_analysis"=1,
  #   "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "1_Group_log_pearson_", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Group",
  #                                  "covariates"="",
  #                                  "family_test"=c("wilcoxon"),
  #                                  "transformation"="quantile_5",
  #                                  "depth_analysis"=1,
  #                                  "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "1_Group_quantile_5_wilcoxon_", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Group",
  #                                  "covariates"="",
  #                                  "family_test"=c("wilcoxon"),
  #                                  "transformation"="quantile_5",
  #                                  "depth_analysis"=3,
  #                                  "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers="DELTAS", areas="PROBE")
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Group_quantile_5_wilcoxon_", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates" =c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="scale",
  #   "depth_analysis"=1,
  #   "filter_p_value" = FALSE)
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers="DELTAS", areas="GENE")
  #
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "1_Phenotest_scale_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="scale",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  #
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, figures="HYPER", markers="DELTAS", areas="GENE")
  #
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_scale_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="log10",
  #   "depth_analysis"=2,
  #   "filter_p_value" = FALSE)
  #
  #
  # # inference_details,result_folder, maxResources, parallel_strategy
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas="GENE")
  #
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "2_Phenotest_log10_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)
  #
  # ####################################################################################
  #
  # inference_details <- expand.grid("independent_variable"= "Phenotest",
  #   "covariates"=c("Covariates1+Covariates2"),
  #   "family_test"=c("gaussian"),
  #   "transformation"="logn",
  #   "depth_analysis"=3,
  #   "filter_p_value" = FALSE)
  #
  # areas_selection <- rownames(signal_data)[1:100]
  # # inference_details,result_folder, maxResources, parallel_strategy
  # semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy, areas_selection=areas_selection, areas="PROBE")
  #
  # fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_logn_gaussian_Covariates1_Covariates2", extension = "csv")
  # localFileRes <- read.table(fileToRead, sep=";")
  # testthat::expect_true(nrow(localFileRes)>0)


  ####################################################################################
  unlink(tempFolder,recursive = TRUE)
  semseeker:::close_env()

  })

