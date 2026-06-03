analyze_batch <- function(signal_data, sample_sheet)
{

  ssEnv <- get_session_info()
  batch_id <- ssEnv$running_batch_id
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " working on batch:", batch_id, " of ", nrow(signal_data), " rows and ", ncol(signal_data), " samples.")
  colnames(signal_data) <- name_cleaning(colnames(signal_data))
  # Keep Sample_ID in sync with name_cleaning so sample_group_check passes
  sample_sheet$Sample_ID <- name_cleaning(sample_sheet$Sample_ID)
  # AI-027: read via unified dispatcher; CASE 2 (streaming merge) lets
  # the SEM step pick up raw bed/bedgraph files when the SIGNAL_MEAN
  # pivot has not been materialised yet.
  signal_pivot <- read_pivot("SIGNAL", "MEAN", "POSITION", "WHOLE")
  if (is.null(signal_pivot)) {
    # Transparent conversion: WGBS/LONGREAD coordinate input → synthetic probe IDs
    signal_data <- normalize_signal_input(signal_data)
    signal_data <- substitute_infinite(signal_data)
    signal_data <- inpute_missing_values(signal_data)
  } else
  {
    signal_data <- as.data.frame(signal_pivot$collect())
    if("CHR" %in% colnames(signal_data))
      signal_data <- position_pivot_to_probe(signal_data)
  }

  signal_data <- as.data.frame(signal_data)
  ssEnv <- get_meth_tech(signal_data)
  ssEnv <- get_session_info()

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " I will work on:", nrow(signal_data), " PROBES.")

  # Build probe features — path depends on technology:
  #   WGBS / LONGREAD: coordinates are encoded in the synthetic probe IDs;
  #                    no Bioconductor annotation package is needed.
  #   Illumina (K850/K450/K27): use Bioconductor annotation as before.
  if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {
    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " building probe_features from synthetic probe IDs (", ssEnv$tech, ")")
    probe_features <- coord_probe_features(rownames(signal_data))
    probe_features <- sort_by_chr_and_start(probe_features)
    signal_data    <- signal_data[match(probe_features$PROBE, rownames(signal_data)), ]
  } else {
    probe_features <- probe_features_get("PROBE")
    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
              " loaded probe_features from Bioconductor annotation")
    probe_features <- probe_features[(probe_features$PROBE %in% rownames(signal_data)), ]
    probe_features <- sort_by_chr_and_start(probe_features)
    signal_data    <- signal_data[rownames(signal_data) %in% probe_features$PROBE, ]
    probe_features <- probe_features[probe_features$PROBE %in% rownames(signal_data), ]
    signal_data    <- signal_data[match(probe_features$PROBE, rownames(signal_data)), ]
  }

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
  log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"), " post-signal_save  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

  # Reference population subset: probe IDs are kept as rownames (preserved by
  # data.frame column subsetting in R). NO extra "PROBE" column wrapper — the
  # downstream signal_range_values() reads probe IDs from rownames.
  referencePopulationSampleSheet <- sample_sheet[sample_sheet$Sample_Group == "Reference", ]
  log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"), " pre-subset       mem_MB=", round(sum(gc()[, "(Mb)"]), 1), " n_ref=", nrow(referencePopulationSampleSheet))
  referencePopulationMatrix <- signal_data[, referencePopulationSampleSheet$Sample_ID, drop = FALSE]
  log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"), " post-subset      mem_MB=", round(sum(gc()[, "(Mb)"]), 1), " dim=", nrow(referencePopulationMatrix), "x", ncol(referencePopulationMatrix))

  if (plyr::empty(referencePopulationMatrix) || ncol(referencePopulationMatrix) < 1) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Empty signal_data ", format(Sys.time(), "%a %b %d %X %Y"))
    stop()
  }

  log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"), " pre-thresholds   mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
  populationControlRangeBetaValues <- as.data.frame(signal_range_values(referencePopulationMatrix,batch_id, probe_features))
  log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"), " post-thresholds  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
  rm(referencePopulationMatrix)
  gc()
  log_event("DEBUG_MEM: ", format(Sys.time(), "%a %b %d %X %Y"), " post-rm-refmatr  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

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
    # Use intersect to read column names WITHOUT allocating a subset of signal_data.
    # The previous form `colnames(signal_data[, populationSampleSheet$Sample_ID])`
    # allocated a full data.frame copy (~6 GB on ewas_data_hub Case/Control) only
    # to read its names — discarded immediately. Pure memory waste under big inputs.
    populationMatrixColumns <- intersect(colnames(signal_data), populationSampleSheet$Sample_ID)

    if (length(populationMatrixColumns)==0) {
      log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Population ",sample_group, " is empty, probably the samples of this group are present in another group ? ")
    }
    else
    {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Working on population ",sample_group, " with ", nrow(signal_data), " probes.")
      # Subset once into a named temp variable so it can be freed explicitly after the call.
      population_signal <- signal_data[, populationMatrixColumns, drop = FALSE]
      analyze_population(
        signal_data = population_signal,
        sample_sheet = populationSampleSheet,
        signal_thresholds = populationControlRangeBetaValues,
        probe_features = probe_features
      )
      rm(population_signal)
      gc()
    }
  }
  if (exists("signal_data", envir = environment(), inherits = FALSE))
    rm("signal_data", envir = environment())
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Batch completed:", batch_id)

}
