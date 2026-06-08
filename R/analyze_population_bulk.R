#' analyze_population_bulk — vectorized population analysis (AI-042, 2026-06-08)
#'
#' Drop-in replacement per analyze_population() che SOSTITUISCE il per-sample
#' loop (con dump bedgraph per ogni sample x marker x figure) con operazioni
#' bulk Polars sul SIGNAL pivot wide gia' scritto da signal_save.
#'
#' Pipeline:
#'   1. SIGNAL pivot esiste gia' (signal_save) -> niente per-sample SIGNAL dump.
#'   2. DELTAS_HYPER pivot = (SIGNAL - high_threshold) clip lower=0  (per-col bulk)
#'      DELTAS_HYPO  pivot = (low_threshold - SIGNAL) clip lower=0
#'   3. MUTATIONS_{HYPER,HYPO} pivot = (DELTAS > 0) cast Int32
#'   4. DELTAR_{HYPER,HYPO}    pivot = DELTAS / (high - low)
#'   5. LESIONS_{HYPER,HYPO}   pivot = sliding-window enrichment binomial test
#'      per sample (operando sulle colonne del MUTATIONS pivot, in-memory).
#'
#' Output identici alle versioni per-sample bit-per-bit per MUTATIONS/LESIONS
#' (rank-invariante sotto monotona, vedi memory option-A-beta-to-m-upstream)
#' e numericamente equivalenti per DELTAS/DELTAR (stessa formula, applicata
#' bulk invece di per sample).
#'
#' Risparmio atteso su ewas (4013 sample):
#'   - per-sample loop: ~11 h (write+gzip bedgraph + decompress + read.delim)
#'   - bulk: ~10-30 min (lazy ops Polars + per-sample LESIONS in memoria)
#'
#' @param signal_data passato per compatibilita' di firma, ignorato (usiamo
#'   il pivot SIGNAL gia' su disco da signal_save).
#' @param sample_sheet df con Sample_ID, Sample_Group; usato per liste sample.
#' @param signal_thresholds df con CHR, START, END,
#'   signal_superior_thresholds (high), signal_inferior_thresholds (low),
#'   signal_median_values.
#' @param probe_features passato per compatibilita'.
#'
#' @keywords internal
#' @noRd
analyze_population_bulk <- function(signal_data, sample_sheet,
                                    signal_thresholds, probe_features) {

  ssEnv <- get_session_info()
  start_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [analyze_population_bulk] start (AI-042 vectorized)")

  # ---- Step 1: prepare paths and inputs ------------------------------------
  pivots_dir <- file.path(ssEnv$result_folderData, "Pivots")
  dir.create(pivots_dir, recursive = TRUE, showWarnings = FALSE)
  for (m in c("SIGNAL", "DELTAS", "DELTAR", "MUTATIONS", "LESIONS")) {
    dir.create(file.path(pivots_dir, m), recursive = TRUE, showWarnings = FALSE)
  }

  signal_pivot_path <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION", "WHOLE")
  if (!file.exists(signal_pivot_path)) {
    stop("[analyze_population_bulk] SIGNAL POSITION pivot mancante: ",
         signal_pivot_path,
         " — signal_save() deve essere chiamato prima.")
  }

  # signal_thresholds e' un data.frame in memoria. Polars 1.11 R ha issue con
  # with_columns(cast) sulle lazy frame da data.frame quando CHR e' factor
  # (resta Categorical anche dopo cast). Workaround: scrivere via arrow a
  # parquet temporaneo (arrow rispetta character -> String), poi scan_parquet.
  thr_parquet_path <- tempfile(pattern = "thresholds_bulk_", fileext = ".parquet")
  on.exit(unlink(thr_parquet_path, force = TRUE), add = TRUE)
  arrow::write_parquet(
    data.frame(
      CHR   = as.character(signal_thresholds$CHR),
      START = as.integer(signal_thresholds$START),
      END   = as.integer(signal_thresholds$END),
      .HIGH = signal_thresholds$signal_superior_thresholds,
      .LOW  = signal_thresholds$signal_inferior_thresholds,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    thr_parquet_path
  )
  thr_lazy <- polars::pl$scan_parquet(thr_parquet_path)

  # Carica SIGNAL pivot lazy + cast CHR to String (signal_save lascia Categorical
  # mentre thr_lazy ha String -> mismatch nel join) + join thresholds
  signal_lazy <- read_pivot("SIGNAL", "MEAN", "POSITION", "WHOLE")$
    with_columns(polars::pl$col("CHR")$cast(polars::pl$String))$
    join(thr_lazy, on = c("CHR", "START", "END"), how = "inner")

  # Identifica le colonne sample (escludendo CHR/START/END/.HIGH/.LOW)
  schema_cols <- names(signal_lazy$collect_schema())
  coord_cols  <- c("CHR", "START", "END", ".HIGH", ".LOW")
  sample_cols <- setdiff(schema_cols, coord_cols)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [analyze_population_bulk] joined SIGNAL with thresholds: ",
            length(sample_cols), " sample columns, ",
            nrow(signal_thresholds), " positions")

  # ---- Step 2: DELTAS_HYPER, DELTAS_HYPO bulk -----------------------------
  # DELTAS_HYPER[s, p] = max(SIGNAL[s,p] - high[p], 0)
  # DELTAS_HYPO[s, p]  = max(low[p] - SIGNAL[s,p], 0)
  for (figure in c("HYPER", "HYPO")) {
    deltas_path <- pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    if (file.exists(deltas_path)) {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] DELTAS_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    deltas_exprs <- lapply(sample_cols, function(s) {
      if (figure == "HYPER") {
        diff_expr <- polars::pl$col(s)$sub(polars::pl$col(".HIGH"))
      } else {
        diff_expr <- polars::pl$col(".LOW")$sub(polars::pl$col(s))
      }
      diff_expr$clip(lower_bound = 0)$alias(s)
    })
    out_lazy <- do.call(signal_lazy$with_columns, deltas_exprs)$
      drop(c(".HIGH", ".LOW"))
    out_lazy$sink_parquet(deltas_path)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] DELTAS_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
  }

  # ---- Step 3: MUTATIONS_{HYPER,HYPO} bulk = (DELTAS > 0) cast Int32 -----
  for (figure in c("HYPER", "HYPO")) {
    mut_path <- pivot_file_name_parquet("MUTATIONS", figure, "POSITION", "WHOLE")
    if (file.exists(mut_path)) {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] MUTATIONS_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    deltas_path <- pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    deltas_lazy <- polars::pl$scan_parquet(deltas_path)
    mut_exprs <- lapply(sample_cols, function(s) {
      polars::pl$col(s)$gt(0)$cast(polars::pl$Int32)$alias(s)
    })
    out_lazy <- do.call(deltas_lazy$with_columns, mut_exprs)
    out_lazy$sink_parquet(mut_path)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] MUTATIONS_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
  }

  # ---- Step 4: DELTAR_{HYPER,HYPO} bulk = DELTAS / (high - low) ----------
  for (figure in c("HYPER", "HYPO")) {
    deltar_path <- pivot_file_name_parquet("DELTAR", figure, "POSITION", "WHOLE")
    if (file.exists(deltar_path)) {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] DELTAR_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    deltas_path <- pivot_file_name_parquet("DELTAS", figure, "POSITION", "WHOLE")
    deltas_lazy <- polars::pl$scan_parquet(deltas_path)$
      with_columns(polars::pl$col("CHR")$cast(polars::pl$String))$
      join(thr_lazy$select(c("CHR", "START", "END", ".HIGH", ".LOW")),
           on = c("CHR", "START", "END"), how = "inner")
    # dividend = high - low (avoid 0 by adding tiny eps a la deltar_single_sample)
    deltar_exprs <- lapply(sample_cols, function(s) {
      div_expr <- polars::pl$col(".HIGH")$sub(polars::pl$col(".LOW"))$
        clip(lower_bound = 1e-9)
      polars::pl$col(s)$truediv(div_expr)$alias(s)
    })
    out_lazy <- do.call(deltas_lazy$with_columns, deltar_exprs)$
      drop(c(".HIGH", ".LOW"))
    out_lazy$sink_parquet(deltar_path)
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] DELTAR_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " sec")
  }

  # ---- Step 5: LESIONS_{HYPER,HYPO} per-sample column processing ---------
  # LESIONS calculation: AI-044 (2026-06-08) refactor a physical bp-window
  # via lesions_get(). Soppianta la vecchia logica row-count (sliding_window_size)
  # commentata sotto per reference; la finestra fisica e' biologicamente piu'
  # rigorosa per il concetto "aggregati mono-direzionali localizzati".
  window_kbp     <- as.numeric(ssEnv$lesion_window_kbp)
  bonf_threshold <- as.numeric(ssEnv$bonferroni_threshold)
  CHUNK_SAMPLES  <- 200L  # gruppi di sample per limitare RAM

  for (figure in c("HYPER", "HYPO")) {
    lesions_path <- pivot_file_name_parquet("LESIONS", figure, "POSITION", "WHOLE")
    if (file.exists(lesions_path)) {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " [bulk] LESIONS_", figure, " pivot already exists, skip")
      next
    }
    t0 <- Sys.time()
    mut_path <- pivot_file_name_parquet("MUTATIONS", figure, "POSITION", "WHOLE")

    sample_chunks <- split(sample_cols,
                          ceiling(seq_along(sample_cols) / CHUNK_SAMPLES))
    tmp_dir <- tempfile("lesions_chunks_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    for (ci in seq_along(sample_chunks)) {
      cols_this <- sample_chunks[[ci]]
      # Carica solo le colonne di questo chunk + coordinate (ordinato per CHR/START)
      mut_chunk <- polars::pl$scan_parquet(mut_path)$
        select(c("CHR", "START", "END", cols_this))$
        sort(c("CHR", "START"))$
        collect()
      mut_df <- as.data.frame(mut_chunk)
      mut_df$CHR <- as.character(mut_df$CHR)

      # Delega LESIONS computation a lesions_new (physical bp window)
      les_mat <- lesions_get(mut_df, cols_this,
                             window_kbp = window_kbp,
                             bonf_threshold = bonf_threshold)

      les_df <- data.frame(
        CHR   = mut_df$CHR,
        START = mut_df$START,
        END   = mut_df$END,
        les_mat,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      chunk_file <- file.path(tmp_dir, sprintf("%03d.parquet", ci))
      polars::as_polars_df(les_df)$write_parquet(chunk_file)
      rm(mut_chunk, mut_df, les_mat, les_df); invisible(gc(verbose = FALSE))
    }

    # Concat dei chunk per colonna (tutti hanno stesso ordine righe).
    chunk_paths <- list.files(tmp_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(chunk_paths) == 1L) {
      file.copy(chunk_paths[1], lesions_path, overwrite = TRUE)
    } else {
      les_combined <- polars::pl$scan_parquet(chunk_paths[1])
      for (i in 2:length(chunk_paths)) {
        les_next <- polars::pl$scan_parquet(chunk_paths[i])
        les_combined <- les_combined$join(les_next,
                                          on = c("CHR", "START", "END"),
                                          how = "inner")
      }
      les_combined$sink_parquet(lesions_path)
    }

    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " [bulk] LESIONS_", figure, " written in ",
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
              " min (",  length(sample_chunks), " chunks of ",
              CHUNK_SAMPLES, " samples, window_kbp=", window_kbp, ")")
  }

  # ============================================================
  # LEGACY ROW-WINDOW LESIONS CODE (AI-042 v1, soppiantata da
  # lesions_new in AI-044 il 2026-06-08). Commentato per reference;
  # NON eseguito. Eliminabile in cleanup futuro.
  # ============================================================
  # Algoritmo legacy (row-count window via sliding_window_size):
  #   ENRICHMENT      = rolling_sum(MUT, sliding_window_size, center)
  #   BASEPAIR_COUNT  = rolling_max(START) - rolling_min(START) (variabile per riga)
  #   p0              = sum(MUT) / length(MUT) per-CHR
  #   lesionpValue    = pbinom(ENRICHMENT-1, size=W, prob=p0, lower=FALSE)
  #   LESIONS_BOOL    = lesionpValue < bonferroni / (n_probes * log10(BASEPAIR_COUNT))
  #
  # Differenza chiave vs lesions_new: window in N righe vs ±X kbp fisici.
  # Il row-window non e' biologicamente coerente perche' span fisico varia
  # 400x fra CpG islands e regioni intergeniche.

  end_time <- Sys.time()
  total_min <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " [analyze_population_bulk] completed in ", total_min, " min")
  invisible(NULL)
}
