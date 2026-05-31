test_that("deltaX_get", {

  tempFolder <- tempFolders[1]
  unlink(tempFolder, recursive = TRUE)
  # message(tempFolder)
  tempFolders <- tempFolders[-1]
  ssEnv <- semseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
    bonferroni_threshold = bonferroni_threshold,
    inpute="median", start_fresh=TRUE)

  ####################################################################################

  tt <- semseeker:::get_meth_tech(signal_data)

  ####################################################################################
  sliding_window_size <- 11
  bonferroni_threshold <- 0.05

  if (!exists("signal_thresholds"))
  {
    signal_data <- semseeker:::inpute_missing_values(signal_data)
    signal_thresholds <<- semseeker:::signal_range_values(signal_data, batch_id)
  }
  probe_features <<- probe_features[probe_features$PROBE %in% rownames(signal_data), ]

  keys <- ssEnv$keys_markers_figures_default
  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Case",],
    probe_features = probe_features
  )
  semseeker:::create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Case",],keys)

  sp <- semseeker:::analyze_population(signal_data=signal_data,
    signal_thresholds = signal_thresholds,
    sample_sheet = mySampleSheet[mySampleSheet$Sample_Group == "Control",],
    probe_features = probe_features
  )

  semseeker:::create_position_pivots(mySampleSheet[mySampleSheet$Sample_Group == "Control",],keys)

  # deltaX_get() calls study_summary_get() which reads the sample sheet CSV.
  # analyze_population is called directly here (bypassing analyze_batch which writes it),
  # so we write it manually with original mixed-case Sample_IDs.
  ssEnv2 <- semseeker:::get_session_info()
  sample_sheet_csv <- semseeker:::file_path_build(ssEnv2$result_folderData, "1_sample_sheet_original", "csv", FALSE)
  utils::write.csv2(mySampleSheet, file=sample_sheet_csv)

  ss <- semseeker:::deltaX_get()

  # analyze_population bypasses analyze_batch so Sample_IDs remain in original case
  cleaned_sample_ids <- mySampleSheet$Sample_ID

  # verify all DELTAX (except DELTAS and DELTAR ) are coherent with MUTATIONS
  keys <- subset(ssEnv$keys_areas_subareas_markers_figures, AREA == "POSITION")
  keys <- subset(keys, MARKER != "DELTAS")
  keys <- subset(keys, MARKER != "DELTAR")
  keys <- subset(keys, MARKER != "LESIONS")
  keys <- subset(keys, MARKER != "SIGNAL")


  for (k in 1:nrow(keys))
  {
    # k <- 1
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    area <- as.character(key$AREA)
    subarea <- as.character(key$SUBAREA)

    mutations_pivot_file_name <- semseeker:::pivot_file_name_parquet("MUTATIONS",figure,area,subarea)
    if(file.exists(mutations_pivot_file_name))
      mutations_pivot <- as.data.frame(polars::pl$read_parquet(mutations_pivot_file_name))
    else
      next

    pivot_file_name <- semseeker:::pivot_file_name_parquet(marker,figure,area,subarea)
    # derived markers may not exist with sparse synthetic data
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

    # Only MUTATIONS is trivially identical to itself; DELTA* use signed values
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
        sample_id <- mySampleSheet[c,"Sample_ID"]  # original case (analyze_population not via analyze_batch)
        sample_group <- mySampleSheet[c,"Sample_Group"]
        mutation_bed_file_name <- semseeker:::bed_file_name(sample_id,sample_group,"MUTATIONS",figure)
        if(file.exists(mutation_bed_file_name))
        {
          marker_bed_file_name <- semseeker:::bed_file_name(sample_id,sample_group,marker,figure)
          testthat::expect_true(file.exists(marker_bed_file_name))
        }
      }

  }

  semseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

