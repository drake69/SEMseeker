test_that("get_result_inference", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]

  ####################################################################################


  semseeker:::semseeker( sample_sheet =  mySampleSheet,signal_data =  signal_data,result_folder = tempFolder,parallel_strategy=parallel_strategy, figures="HYPER",
    markers=c("DELTAQ"), areas=c("PROBE","CHR","GENE"), start_fresh = TRUE, inpute="median")

  ####################################################################################

  inferenceFolder <- file.path(tempFolder,"Inference")
  inference_details <- expand.grid("independent_variable"= "Phenotest",
    "covariates"=c("Covariates1+Covariates2"),
    "family_test"=c("gaussian"),
    "transformation"="",
    "depth_analysis"=3,
    "filter_p_value" = FALSE)
  semseeker:::association_analysis(inference_details = inference_details, result_folder = tempFolder, parallel_strategy=parallel_strategy,figures="HYPER",
    markers=c("DELTAS","DELTAQ"), areas="GENE")
  fileToRead <- semseeker:::file_path_build(inferenceFolder, "3_Phenotest_gaussian_Covariates1_Covariates2", extension = "csv")
  localFileRes <- read.table(fileToRead, sep=";")
  testthat::expect_true(nrow(localFileRes)>0)

})
