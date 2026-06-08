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
  # Per-chromosome sort+sink (era sort+sink globale → Jetsam OOM kill silenzioso
  # su 64GB Mac quando la lazy chain join+sort+collect superava soglie macOS
  # memorystatus). Polars 1.11 R bindings non hanno streaming-mode + projection
  # pushdown attraverso join è incompleta, quindi gestiamo le chr dal lato R
  # (probe_features è già una data.frame piccola in RAM) e iteriamo ognuna come
  # query lazy indipendente filtrata sul subset di PROBE corrispondente.
  pf <- probe_features_get("PROBE")    # data.frame: PROBE, CHR, START, END
  chrs_all <- as.character(pf$CHR)
  chr_order_key <- function(ch) {
    cl <- sub("^chr", "", ch)
    n  <- suppressWarnings(as.integer(cl))
    ifelse(!is.na(n), n,
           ifelse(cl == "X", 23L,
                  ifelse(cl == "Y", 24L,
                         ifelse(cl %in% c("M", "MT"), 25L, 99L))))
  }
  chrs <- unique(chrs_all)
  chrs <- chrs[order(chr_order_key(chrs))]
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"),
            " pre-chunked-sort mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
            " (chromosomes: ", length(chrs), ")")

  tmp_dir <- tempfile("signal_chr_chunks_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  pp_lazy <- polars::as_polars_df(pf)$lazy()
  # Scan parquet UNA SOLA VOLTA fuori dal loop: la lazy frame è immutabile,
  # ogni $filter/$join in iterazione produce una nuova lazy frame senza
  # ri-aprire il parquet. Riduce overhead per iter su big matrix.
  sd_lazy <- read_pivot("SIGNAL", "MEAN", "PROBE", "WHOLE")$
              with_columns(polars::pl$col("AREA")$alias("PROBE"))

  for (i in seq_along(chrs)) {
    ch <- chrs[i]
    chunk_file <- file.path(tmp_dir, sprintf("%03d_%s.parquet", i, ch))
    # Lazy chain per chunk: pp filtered → join sd_lazy → drop → sort → sink.
    # Il pushdown del filter `pp_lazy$filter(CHR == ch)` (~25k sonde) attraverso
    # l'inner join restringe sd_lazy alle sole righe di quel chr.
    # NON usiamo `is_in(probes_chr)` su sd_chr: in polars R 1.11 il vettore R
    # viene interpretato come nomi di colonna invece che valori → "not found".
    sd_chr <- pp_lazy$filter(polars::pl$col("CHR") == ch)$
                      join(sd_lazy, on = "PROBE", how = "inner")$
                      drop(c("PROBE", "PROBE_WHOLE", "AREA"))$
                      sort(c("START", "END"), descending = FALSE)
    sd_chr$sink_parquet(chunk_file)
    # Defensive cleanup: rilascia R-side reference, forza gc() per evitare
    # accumulo di stato Polars (mmap, cache, plan) fra iterazioni su big matrix.
    rm(sd_chr); invisible(gc(verbose = FALSE))
  }
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"),
            " post-chunked-sort mem_MB=", round(sum(gc()[, "(Mb)"]), 1),
            " (wrote ", length(chrs), " chunks)")

  # Libera pp_lazy + sd_lazy + pf prima del concat — la lazy frame del concat
  # streamerà i chunk dal filesystem, non serve mantenere queste references.
  rm(sd_lazy, pp_lazy, pf); invisible(gc(verbose = FALSE))

  chunk_paths <- list.files(tmp_dir, pattern = "\\.parquet$", full.names = TRUE)
  # pl$concat() vuole `...` varargs di LazyFrame, NON una lista, quindi do.call
  # per unpack-are. Output finale = concat di chunk già pre-ordinati per chr
  # (ordine canonico via chr_order_key) → equivalente al sort globale CHR/START/END.
  do.call(polars::pl$concat,
          lapply(chunk_paths, polars::pl$scan_parquet))$
    sink_parquet(pivot_file_name_pos)
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-sink-position mem_MB=", round(sum(gc()[, "(Mb)"]), 1))

  rm(signal_data, pp)
  gc()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), "Saved signal data.")
  log_event("DEBUG_MEM_SS: ", format(Sys.time(), "%a %b %d %X %Y"), " post-rm-lazyframes   mem_MB=", round(sum(gc()[, "(Mb)"]), 1))
}
