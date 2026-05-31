test_that("sort_by_chr_and_start", {


  tempFolder <- tempFolders[1]
  tempFolders <- tempFolders[-1]
  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")

  #order not matching
  second <- semseeker:::sort_by_chr_and_start( probe_features[order(probe_features$START),])

  testthat::expect_true( semseeker:::test_match_order( semseeker:::sort_by_chr_and_start(probe_features)$ABSOLUTE,second$ABSOLUTE))

  ####################################################################################

  # semseeker:::close_env()
}
)
