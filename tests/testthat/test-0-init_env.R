test_that("test_init_env", {

  # run this test as first other wise the reuse of session is not testable
  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]

  # test minimum number of parameters
  ssEnv <- semseeker:::init_env(tempFolder)
  testthat::expect_true(ssEnv$result_folderData == normalizePath(file.path(tempFolder, "Data"), mustWork = FALSE))
  semseeker:::close_env()

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("PROBE"), subareas= c(""), showprogress=TRUE, start_fresh=TRUE, maxResources = 10)
  # POSITION is always injected internally by keys_create() (needed for study_summary_total).
  # So requesting areas=c("PROBE") yields keys containing PROBE + POSITION only.
  testthat::expect_true(nrow(subset(ssEnv$keys_areas_subareas_markers_figures, !AREA %in% c("PROBE", "POSITION")))==0)


  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas = c("GENE"), subareas= c("BODY"), showprogress=TRUE, start_fresh=TRUE, maxResources = 10)
  testthat::expect_true(ssEnv$maxResources == "10")
  semseeker:::close_env()

  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,areas = c("GENE","DMR"), showprogress=TRUE, start_fresh=TRUE)
  testthat::expect_true(length(unique(ssEnv$keys_areas_subareas[,"AREA"]))!=1)

  unlink(tempFolder,recursive = TRUE)
  assign("ssEnv", NULL, envir=.pkgglobalenv)

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, figures="HYPPO"), "ERROR:")
  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, markers="HYPPO"), "ERROR:")

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, areas="HYPPO"), "ERROR:")

  ####################################################################################

  # expect_error( semseeker:::init_env(tempFolder, parallel_strategy ="cluster"), "ERROR:")

  ####################################################################################

  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, subareas="HYPPO"), "ERROR:")


  testthat::expect_error( semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, test_it="HYPPO"), "ERROR:")

  ####################################################################################

  unlink(tempFolder,recursive = TRUE)
  # semseeker:::close_env()
}
)
