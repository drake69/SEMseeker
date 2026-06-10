#--- signal_range_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of signal values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation (probes × samples)
#' @param batch_id character string identifying the batch; used to name the cached parquet output file
#' @param probe_features data.frame of probe annotations (from probe_features_get) with columns PROBE, CHR, START, END
#'
#' @return data.frame of per-probe thresholds with columns signal_inferior_thresholds,
#'   signal_superior_thresholds, signal_median_values, iqr, q1, q3, PROBE, CHR, START, END
#' @importFrom doRNG %dorng%

signal_range_values <- function(populationMatrix, batch_id, probe_features) {


  ssEnv <- get_session_info()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting signal thresholds calculation.")
  thresholds_file_name <- file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds"),"parquet")
  if(file.exists(thresholds_file_name))
  {
    result <- as.data.frame(polars::pl$read_parquet(thresholds_file_name))
    return(result)
  }
  if (sum(is.na(populationMatrix)) > 0)
  {
    msg <- paste0("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), " There are missing values in the population matrix, apply the parameter inpute or remove the missing values.")
    log_event(msg)
    stop(msg)
  }

  # Vectorized threshold computation via matrixStats (C-level, ~10x faster
  # than row-wise apply). Scales to 28M positions (whole-genome WGBS).
  # Caller passes a probe-indexed data.frame or matrix (no PROBE column);
  # probe IDs live in rownames. as.matrix() forces numeric matrix needed by
  # matrixStats but is a no-op if input is already a matrix.
  mat <- if (is.matrix(populationMatrix)) populationMatrix else as.matrix(populationMatrix)
  iqr_times <- as.numeric(ssEnv$iqrTimes)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Computing thresholds for ", nrow(mat), " positions x ",
            ncol(mat), " samples (vectorized).")

  q1  <- matrixStats::rowQuantiles(mat, probs = 0.25, na.rm = TRUE)
  med <- matrixStats::rowMedians(mat, na.rm = TRUE)
  q3  <- matrixStats::rowQuantiles(mat, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  rmins <- matrixStats::rowMins(mat, na.rm = TRUE)
  rmaxs <- matrixStats::rowMaxs(mat, na.rm = TRUE)

  # Thresholds clamped to observed min/max
  signal_inferior <- pmax(q1 - iqr_times * iqr, rmins)
  signal_superior <- pmin(q3 + iqr_times * iqr, rmaxs)

  result <- data.frame(
    signal_inferior_thresholds = signal_inferior,
    signal_superior_thresholds = signal_superior,
    signal_median_values       = med,
    iqr                        = iqr,
    q1                         = q1,
    q3                         = q3
  )
  rm(mat, q1, med, q3, iqr, rmins, rmaxs, signal_inferior, signal_superior)
  gc()

  #
  colnames(result) <- c("signal_inferior_thresholds","signal_superior_thresholds", "signal_median_values","iqr","q1","q3")
  result$PROBE <- rownames(populationMatrix)
  if(nrow(result) != nrow(populationMatrix))
    stop("I'M STOPPING HERE, No thresholds defined for the population.")

  result <- polars::as_polars_df(result)$lazy()
  probe_features <- polars::as_polars_df(probe_features)$lazy()
  # rename AREA as PROBE
  result <- probe_features$join(
    result,
    on = c("PROBE"),
    how = "inner"
  )
  # Drop PROBE_WHOLE before writing thresholds: it's the canonical
  # (AREA, SUBAREA) name for the probe-level baseline (kept on
  # probe_features so downstream association code can do dynamic
  # `probe_features[[area_subarea]]` lookup) but it's a 1:1 alias of
  # PROBE, so persisting it to the thresholds parquet is pure waste.
  if ("PROBE_WHOLE" %in% names(result))
    result <- result$drop(c("PROBE_WHOLE"))
  result <- result$sort(c("CHR", "START","END"), descending = FALSE)
  result <- result$collect()
  result$write_parquet(thresholds_file_name)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Thresholds defined for: ", nrow(result), " probe_features.")
  gc()
  return(as.data.frame(result))
}
