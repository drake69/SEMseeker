signal_save <- function(signal_data, sample_sheet, batch_id)
{
  ssEnv <- get_session_info()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saving signal data.")
  pivot_file_name_pos <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION", "WHOLE")
  if (file.exists(pivot_file_name_pos)) {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Signal data already saved.")
    return()
  }

  signal_data <- signal_data[, unique(sample_sheet$Sample_ID), drop = FALSE]

  # ------------------------------------------------------------------
  # WGBS / LONGREAD path — coordinates are encoded in synthetic probe IDs.
  # No Bioconductor annotation join is needed.
  # ------------------------------------------------------------------
  if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              "Saving signal data for ", ssEnv$tech,
              " (extracting coordinates from synthetic probe IDs).")

    # Save probe-indexed parquet (synthetic IDs as AREA column)
    signal_probe           <- signal_data
    signal_probe$AREA      <- rownames(signal_data)
    signal_probe           <- signal_probe[, c(ncol(signal_probe),
                                               seq_len(ncol(signal_probe) - 1))]
    pivot_file_name_probe  <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")
    polars::as_polars_df(signal_probe)$write_parquet(pivot_file_name_probe)
    rm(signal_probe)

    # Build position-indexed parquet directly from synthetic probe IDs
    coords       <- probe_id_to_coord(rownames(signal_data))
    signal_pos   <- data.frame(
      CHR   = coords$CHR,
      START = coords$START,
      END   = coords$END,
      signal_data,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    # Sort by chromosome (natural order) then position
    chr_order    <- order(
      suppressWarnings(as.integer(signal_pos$CHR)),   # numeric chrs first
      nchar(signal_pos$CHR),                          # X/Y/M after
      signal_pos$CHR,
      signal_pos$START
    )
    signal_pos   <- signal_pos[chr_order, ]
    polars::as_polars_df(signal_pos)$write_parquet(pivot_file_name_pos)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data (", ssEnv$tech, ").")
    gc()
    return()
  }

  # ------------------------------------------------------------------
  # Illumina path — join with Bioconductor annotation to get CHR/START/END
  # ------------------------------------------------------------------
  signal_data$AREA <- rownames(signal_data)
  signal_data      <- signal_data[, c(ncol(signal_data), seq_len(ncol(signal_data) - 1))]

  pivot_file_name_probe <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")
  polars::as_polars_df(signal_data)$write_parquet(pivot_file_name_probe)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Signal data saved with probe.")

  rm(signal_data)
  gc()
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-probe-write+gc  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

  # AI-027: read via unified dispatcher. The PROBE pivot was just
  # written above (line 62), so CASE 1 (cached parquet) is always taken.
  signal_data <- read_pivot("SIGNAL", "MEAN", "PROBE", "WHOLE")
  pp          <- polars::as_polars_df(probe_features_get("PROBE"))$lazy()
  signal_data <- signal_data$with_columns(polars::pl$col("AREA")$alias("PROBE"))
  signal_data <- pp$join(signal_data, on = "PROBE", how = "inner")
  signal_data <- signal_data$drop(c("PROBE", "PROBE_WHOLE", "AREA"))
  signal_data <- signal_data$sort(c("CHR", "START", "END"), descending = FALSE)
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " pre-collect-position mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
  signal_data$collect()$write_parquet(pivot_file_name_pos)
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-position-write  mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

  rm(signal_data, pp)
  gc()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data.")
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-rm-lazyframes   mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
}
