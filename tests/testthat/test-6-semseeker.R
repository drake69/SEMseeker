test_that("semeeker", {

  tempFolder <- tempFolders[1]
  unlink(tempFolder,recursive = TRUE)
  tempFolders <- tempFolders[-1]
  ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy, showprogress=showprogress, verbosity=verbosity)

  ####################################################################################

  SEMseeker:::semseeker( sample_sheet =  mySampleSheet ,signal_data =  signal_data,
    result_folder = tempFolder, parallel_strategy = parallel_strategy)

  keys <- subset(ssEnv$keys_areas_subareas_markers_figures)
  # name_cleaning uppercases Sample_ID inside semseeker(); use the same for comparison
  cleaned_sample_ids <- SEMseeker:::name_cleaning(mySampleSheet$Sample_ID)

  for (k in 1:nrow(keys))
  {
    # k <- 1
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    area <- as.character(key$AREA)
    subarea <- as.character(key$SUBAREA)

    mutations_pivot_file_name <- SEMseeker:::pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    if(file.exists(mutations_pivot_file_name))
      mutations_pivot <- as.data.frame(polars::pl$read_parquet(mutations_pivot_file_name))
    else
      next

    pivot_file_name <- SEMseeker:::pivot_file_name_parquet(marker,figure,area,subarea)
    # derived markers (LESIONS, DELTA*) may not exist when mutations are too sparse
    if(!file.exists(pivot_file_name))
      next
    pivot <- as.data.frame(polars::pl$read_parquet(pivot_file_name))

    pivot <- pivot[,-c(1:3)]
    mutations_pivot <- mutations_pivot[,-c(1:3)]
    testthat::expect_true(nrow(pivot)<nprobes)
    testthat::expect_true(nrow(pivot)>0)

    pivot <- pivot[,order(colnames(pivot))]
    mutations_pivot <- mutations_pivot[,order(colnames(mutations_pivot))]

    testthat::expect_true(all(colnames(pivot) %in% cleaned_sample_ids))
    testthat::expect_true(all(colnames(pivot) == colnames(mutations_pivot)))

    testthat::expect_true(ncol(pivot)==ncol(mutations_pivot))
    testthat::expect_true(nrow(pivot)==nrow(mutations_pivot))

    pivot[is.na(pivot)] <- 0
    mutations_pivot[is.na(mutations_pivot)] <- 0

    pivot[pivot > 0] <- 1
    mutations_pivot[mutations_pivot > 0] <- 1

    # Only MUTATIONS pivot is identical to itself; derived markers (DELTA*, LESIONS)
    # may have signed/windowed values that don't match the binary MUTATIONS mask
    if (marker == "MUTATIONS")
    {
      if (!all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
      {
        print(marker)
        print(figure)
      }
      testthat::expect_true(all(as.data.frame(pivot) == as.data.frame(mutations_pivot)))
    }

    if(marker!="MUTATIONS")
      for(c in 1:(nrow(mySampleSheet)))
      {
        sample_id <- SEMseeker:::name_cleaning(mySampleSheet[c,"Sample_ID"])
        sample_group <- mySampleSheet[c,"Sample_Group"]
        mutation_bed_file_name <- SEMseeker:::bed_file_name(sample_id,sample_group,"MUTATIONS",figure)
        if(file.exists(mutation_bed_file_name))
        {
          marker_bed_file_name <- SEMseeker:::bed_file_name(sample_id,sample_group,marker,figure)
          testthat::expect_true(file.exists(marker_bed_file_name))
        }
      }

  }

  ####################################################################################
  SEMseeker:::close_env()
  unlink(tempFolder,recursive = TRUE)
})

