test_that("set-env-variable", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]

  # sex_chromosome_remove is FALSE by default
  ssEnv <- SEMseeker:::init_env(result_folder= tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), start_fresh=TRUE)
  testthat::expect_false(ssEnv$sex_chromosome_remove)

  # create list arguments
  arguments <- list(areas = c("GENE"), subareas = c("BODY"), sex_chromosome_remove = TRUE)
  arguments  <- SEMseeker:::set_env_variable(arguments =  arguments, var_name =  "sex_chromosome_remove",var_default_value =  TRUE)
  ssEnv  <- SEMseeker:::get_session_info(tempFolder)
  testthat::expect_true(ssEnv$sex_chromosome_remove)
  # check sex_chromosome_remove is not in the list
  testthat::expect_true(is.null(arguments[["sex_chromosome_remove"]]))

  SEMseeker:::close_env()

  # check reopening not fresh is FALSE
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), start_fresh=FALSE)
  testthat::expect_true(ssEnv$sex_chromosome_remove)
  SEMseeker:::close_env()

  # change sex_chromosome_remove to TRUE
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), sex_chromosome_remove=TRUE, start_fresh=FALSE)
  ssEnv  <- SEMseeker:::get_session_info(tempFolder)
  arguments  <- SEMseeker:::set_env_variable(arguments =  arguments, var_name =  "sex_chromosome_remove",var_default_value =  TRUE)
  testthat::expect_true(is.null(arguments[["sex_chromosome_remove"]]))
  ssEnv  <- SEMseeker:::get_session_info(tempFolder)
  testthat::expect_true(ssEnv$sex_chromosome_remove=="TRUE")

  SEMseeker:::close_env()

})
