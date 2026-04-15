test_that("apply_stat_model_sig_formula", {

  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy)

  family_test = "gaussian"
  burdenValue = "dependent_variable"
  independent_variable = "independent_variable"
  covariates = c("cov1","cov2","cov3")

  ff1 <- SEMseeker:::apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
  testthat::expect_true(ff1==as.formula("dependent_variable~independent_variable+cov1+cov2+cov3"))

  family_test = "binomial"
  burdenValue = "dependent_variable"
  independent_variable = "independent_variable"
  covariates = c("cov1","cov2","cov3")

  ff1 <- SEMseeker:::apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
  testthat::expect_true(ff1==as.formula("independent_variable~dependent_variable+cov1+cov2+cov3"))

  family_test = "wilcoxon"
  burdenValue = "dependent_variable"
  independent_variable = "independent_variable"
  covariates = c("cov1","cov2","cov3")

  ff1 <- SEMseeker:::apply_stat_model_sig_formula(family_test, burdenValue, independent_variable, covariates)
  testthat::expect_true(ff1==as.formula("dependent_variable~independent_variable"))


  SEMseeker:::close_env()
  unlink(tempFolder, recursive = T)

})
