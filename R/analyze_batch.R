analyze_batch <- function(signal_data, sample_sheet)
{

  ssEnv <- get_session_info()
  batch_id <- ssEnv$running_batch_id
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " working on batch:", batch_id, " of ", nrow(signal_data), " rows and ", ncol(signal_data), " samples.")
  colnames(signal_data) <- name_cleaning(colnames(signal_data))
  # Keep Sample_ID in sync with name_cleaning so sample_group_check passes
  sample_sheet$Sample_ID <- name_cleaning(sample_sheet$Sample_ID)
  pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION","WHOLE")
  if(!file.exists(pivot_file_name))
  {
    signal_data <- substitute_infinite(signal_data)
    signal_data <- inpute_missing_values(signal_data)
  } else
  {
    signal_data <- as.data.frame(polars::pl$read_parquet(pivot_file_name))
    if("CHR" %in% colnames(signal_data))
      signal_data <- position_pivot_to_probe(signal_data)
  }

  signal_data <- as.data.frame(signal_data)
  ssEnv <- get_meth_tech(signal_data)
  ssEnv <- get_session_info()
  # coverage_analysis(rownames(signal_data))

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work on:", nrow(signal_data), " PROBES.")
  probe_features <- probe_features_get("PROBE")
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " loaded probe_features from Bioconductor annotation")

  probe_features <- probe_features[(probe_features$PROBE %in% rownames(signal_data)),]
  # sort probe features by CHR START and END
  probe_features <- sort_by_chr_and_start(probe_features)
  signal_data <- signal_data[rownames(signal_data) %in% probe_features$PROBE, ]
  probe_features <- probe_features[probe_features$PROBE %in% rownames(signal_data), ]
  signal_data <- signal_data[match(probe_features$PROBE, rownames(signal_data)), ]

  # probe_features <- sort_by_chr_and_start(probe_features)
  if (!test_match_order(row.names(signal_data), probe_features$PROBE)) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Wrong order matching Probes and Methylation data!")
    stop()
  }

  sample_group_checkResult <- sample_group_check(sample_sheet, signal_data)
  if(!is.null(sample_group_checkResult))
  {
    stop(sample_group_checkResult)
  }

  signal_save(signal_data, sample_sheet, batch_id)

  # reference population
  referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  referencePopulationMatrix <- data.frame(PROBE = row.names(signal_data), signal_data[, referencePopulationSampleSheet$Sample_ID])

  # signal_data <- data.frame(PROBE = row.names(signal_data), signal_data[ , which(!(colnames(signal_data)%in%referencePopulationSampleSheet$Sample_ID))]  )

  if (plyr::empty(referencePopulationMatrix) ||
      ncol(referencePopulationMatrix) < 2) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Empty signal_data ", format(Sys.time(), "%a %b %d %X %Y"))
    stop()
  }

  if(!ssEnv$signal_intrasample)
  {
    populationControlRangeBetaValues <- as.data.frame(signal_range_values(referencePopulationMatrix,batch_id, probe_features))
    gc()
  }
  else
  {
    populationControlRangeBetaValues <- NULL
  }

  # remove duplicated samples due to the reference population
  referenceSamples <- sample_sheet[sample_sheet$Sample_Group == "Reference",]
  otherSamples <- sample_sheet[sample_sheet$Sample_Group != "Reference",]
  referenceSamples <- referenceSamples[!(referenceSamples$Sample_ID %in% otherSamples$Sample_ID), ]
  sample_sheet <- rbind(otherSamples, referenceSamples)
  i <- 0
  variables_to_export <- c( "ssEnv", "sample_sheet", "signal_data", "analyze_population",
    "populationControlRangeBetaValues", "probe_features")
  # resultSampleSheet <- foreach::foreach(i = seq_along(ssEnv$keys_sample_groups[,1]), .combine = rbind, .export = variables_to_export ) %dorng%
  for (i in seq_along(ssEnv$keys_sample_groups[,1]))
  {
    sample_group <- ssEnv$keys_sample_groups[i,1]
    populationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == sample_group, ]
    populationMatrixColumns <- colnames(signal_data[, populationSampleSheet$Sample_ID])

    if (length(populationMatrixColumns)==0) {
      log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Population ",sample_group, " is empty, probably the samples of this group are present in another group ? ")
    }
    else
    {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Working on population ",sample_group, " with ", nrow(signal_data), " probes.")
      analyze_population(
        signal_data = signal_data[, populationMatrixColumns],
        sample_sheet = populationSampleSheet,
        signal_thresholds = populationControlRangeBetaValues,
        probe_features = probe_features
      )
      gc()
    }
  }
  if(exists("signal_data"))
    rm(signal_data)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Batch completed:", batch_id)

}
