#' Calculate stochastic epi mutations from a methylation dataset as outcome
#' report of pivot
#'
#' @param signal_data whole matrix of data to analyze.
#' @param sample_sheet name of samplesheet's column to use as control population
#' selector followed by selection value,
#' @param probe_features probe_features detail from 27 to EPIC illumina dataset
#' @param signal_thresholds thresholds defined to calculate epimutations
#' @return files into the result folder with pivot table and bedgraph.
#'   A BANNER is logged once per batch before the per-sample loop showing:
#'   input_positions, beta_range_positions, covered_by_inner_join — allows
#'   immediate audit of cross-run coverage (e.g. Nanopore sample vs Illumina reference).
#' @importFrom doRNG %dorng%
#'
analyze_population <- function(signal_data, sample_sheet,signal_thresholds, probe_features) {

  ssEnv <- get_session_info()
  # #
  start_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " AnalyzePopulation warmingUP ")

  nrow_before <- nrow(signal_data)
  signal_data <- stats::na.omit(signal_data)
  nrow_after <- nrow(signal_data)
  if (nrow_before != nrow_after) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Removed ", nrow_before - nrow_after, " rows with NA values")
    stop()
  }

  ### get signal_values ########################################################
  sample_sheet <- sample_sheet[order(sample_sheet[, "Sample_ID"], decreasing = FALSE), ]
  existent_samples <- colnames(signal_data)
  sample_names <- sample_sheet$Sample_ID
  missed_samples <- setdiff(setdiff(sample_names, existent_samples), "PROBE")

  if (length(missed_samples) != 0) {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " These samples data are missed: ", paste0(missed_samples, sep = " "))
  }

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " WarmedUP AnalyzePopulation")
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Start population analysis")

  # progress_bar <- progress::progress_bar$new(
  #   format = paste("INFO: Performing population analysis [:bar] :percent eta: :eta"),
  #   total = nrow(sample_sheet),
  #   clear = FALSE,
  #   width= 60)

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(sample_sheet))

  variables_to_export <- c("sample_sheet", "signal_data", "analyze_single_sample", "ssEnv",
    "signal_superior_thresholds","deltar_single_sample","signal_inferior_thresholds","iqr","signal_median_values",
    "bt","bonferroni_threshold", "probe_features", "analyze_single_sample_both", "delta_single_sample", "progress_bar",
    "progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","signal_single_sample",
    "get_session_info","bed_file_name","signal_thresholds","update_session_info","normalize_chr")
  i <- 1

  for(i in 1:nrow(sample_sheet)) {
  # foreach::foreach(i =1:nrow(sample_sheet), .export = variables_to_export) %dorng% {
    local_sample_detail <- sample_sheet[i,]
    ssEnv <- get_session_info()
    signal_values <- signal_data[,local_sample_detail$Sample_ID]
    bed_filename <- bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "SIGNAL","MEAN")
    if(!file.exists(bed_filename))
      signal_single_sample( signal_values,local_sample_detail,probe_features)
    if(ssEnv$showprogress)
      progress_bar(sprintf("Saving signal of sample: %s",local_sample_detail$Sample_ID))
  }
  gc()

  # ── Coverage banner — emitted ONCE before per-sample analysis ───────────────
  # Shows how many positions in the current run are covered by signal_thresholds.
  # Especially important for cross-run analysis (e.g. Nanopore sample vs an
  # Illumina reference batch passed via populationControlRangeBetaValues):
  # the user can immediately see if there is a partial or zero overlap.
  # Uses the first sample's signal bed file (just written above) for positions.
  {
    coverage_first_sample <- sample_sheet[1L, ]
    coverage_bed <- bed_file_name(coverage_first_sample$Sample_ID,
                                  coverage_first_sample$Sample_Group,
                                  "SIGNAL", "MEAN")
    if (file.exists(coverage_bed)) {
      coverage_sig <- utils::read.delim(coverage_bed, header = FALSE, sep = "\t")
      colnames(coverage_sig) <- c("CHR", "START", "END", "VALUE")
      coverage_sig$CHR <- normalize_chr(coverage_sig$CHR, "internal")
      coverage_n_input  <- nrow(coverage_sig)
      coverage_n_ranges <- nrow(signal_thresholds)
      coverage_sig_lf <- polars::as_polars_df(
        coverage_sig[, c("CHR", "START", "END")])$
        with_columns(polars::pl$col("CHR")$cast(polars::pl$String))$lazy()
      coverage_thr_lf <- polars::as_polars_df(
        signal_thresholds[, c("CHR", "START", "END")])$
        with_columns(polars::pl$col("CHR")$cast(polars::pl$String))$lazy()
      coverage_n_covered <- coverage_sig_lf$join(
        coverage_thr_lf, on = c("CHR", "START", "END"), how = "inner"
      )$collect()$height
      log_event(
        "BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
        " [analyze_population] Coverage —",
        " input_positions=", coverage_n_input,
        " | beta_range_positions=", coverage_n_ranges,
        " | covered_by_inner_join=", coverage_n_covered,
        " | analysis will run on ", coverage_n_covered, "/", coverage_n_input, " positions"
      )
      rm(coverage_sig, coverage_sig_lf, coverage_thr_lf,
         coverage_first_sample, coverage_bed,
         coverage_n_input, coverage_n_ranges, coverage_n_covered)
    }
  }

  progress_bar <- NULL
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nrow(sample_sheet))

  rm(signal_data)
  # for(i in 1:nrow(sample_sheet)) {
  # .packages loads SEMseeker in each worker so SEMseeker::: lookups resolve.
  # Internal helpers are prefixed with SEMseeker::: because they live in the
  # namespace (not in the caller's frame) and .export does not cover them.
  foreach::foreach(
    i = 1:nrow(sample_sheet),
    .export = variables_to_export,
    .packages = "SEMseeker"
  ) %dorng% {
    # CRITICAL (E-14): multisession workers are fresh R processes where
    # .pkgglobalenv$ssEnv is empty. All internal helpers (bed_file_name,
    # analyze_single_sample, etc.) call get_session_info() which reads from
    # .pkgglobalenv — NOT from the exported `ssEnv` variable. Without this
    # call, multisession workers fail with "get_session_info called without
    # result folder". See engineering-decisions.md §1.3.
    # AI-041: in-memory only; saveRDS would happen N_samples × N_workers
    # times per SEM step otherwise (15 MB per write → catastrophic I/O).
    SEMseeker:::update_session_info(ssEnv, save_to_disk = FALSE)

    local_sample_detail <- sample_sheet[i,]
    bed_filename <- SEMseeker:::bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "SIGNAL","MEAN")
    signal_values <- utils::read.delim(bed_filename, header = FALSE, sep = "\t")
    colnames(signal_values) <- c("CHR", "START", "END", "VALUE")
    signal_values$CHR <- SEMseeker:::normalize_chr(signal_values$CHR, "internal")

    bed_filename <- SEMseeker:::bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "MUTATIONS","HYPER")
    if(!file.exists(bed_filename))
      SEMseeker:::analyze_single_sample( values = signal_values,thresholds = signal_thresholds, figure="HYPER", sample_detail = local_sample_detail)

    bed_filename <- SEMseeker:::bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "MUTATIONS","HYPO")
    if(!file.exists(bed_filename))
      SEMseeker:::analyze_single_sample( values = signal_values,thresholds = signal_thresholds, figure="HYPO", sample_detail = local_sample_detail)

    bed_filename <- SEMseeker:::bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "DELTAS","HYPO")
    if(!file.exists(bed_filename))
      SEMseeker:::delta_single_sample( values = signal_values,thresholds = signal_thresholds , sample_detail = local_sample_detail)

    bed_filename <- SEMseeker:::bed_file_name(local_sample_detail$Sample_ID,local_sample_detail$Sample_Group, "DELTAR","HYPO")
    if(!file.exists(bed_filename))
      SEMseeker:::deltar_single_sample ( values = signal_values, thresholds = signal_thresholds,sample_detail = local_sample_detail)

    if(ssEnv$showprogress)
      progress_bar(sprintf("Performed sample: %s",local_sample_detail$Sample_ID))
  }

  gc()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Row count result:", nrow(sample_sheet))
  if (exists("signal_data", envir = environment(), inherits = FALSE))
    rm("signal_data", envir = environment())

  # AI-041: end-of-batch disk snapshot (workers used save_to_disk=FALSE
  # inside the foreach; here we persist the session exactly once).
  update_session_info(ssEnv, save_to_disk = TRUE)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population analysis ")
  end_time <- Sys.time()
  time_taken <- difftime(end_time,start_time, units = "mins")
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Completed population with summary - Time taken: ", time_taken, " minutes.")

}


