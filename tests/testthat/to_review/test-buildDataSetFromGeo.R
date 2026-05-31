# test_that("build_data_set_from_geo", {
#
#   skip_on_ci()
#   skip_on_covr()
#   skip_on_travis()
#
#
#   sample_sheet <- semseeker::build_data_set_from_geo("GSE132616",tempFolder, 1)
#
#   testthat::expect_true(nrow(mySampleSheet)>0)
#   testthat::expect_true(file.exists(file.path(tempFolder,"/203048410118_R03C01_Red.idat")))
# }
# )
