#' Run depth>1 (area / region level) association for one marker
#'
#' Extracted from association_analysis() (was inline at lines 294-418).
#' Iterates over the area/subarea keys of one marker, reads the
#' corresponding pivot parquet, chunks it (6e6 cells / ncol), transposes
#' and merges with sample_names, then applies the stat model per chunk.
#'
#' @param prep list returned by prepare_study_for_analysis().
#' @param marker character. The marker name (e.g. "MUTATIONS").
#' @param family_test character.
#' @param fileNameResults character. Path of the output CSV.
#' @param filter_p_value logical.
#' @param ssEnv list. Session environment.
#' @param selected_areas character vector or empty.
#' @param results data.frame. Accumulator carried over from depth=1.
#' @param start_time POSIXct. Job start, used by association_analysis_log().
#' @param processed_items integer. Counter carried over from depth=1.
#' @param ... forwarded to apply_stat_model().
#' @return list(results = data.frame, processed_items = integer).
#'   Side effect: writes the CSV via association_analysis_save_results().
#' @keywords internal
run_depth_n_marker <- function(prep, marker, family_test, fileNameResults,
                                filter_p_value, ssEnv, selected_areas,
                                results, start_time, processed_items, ...) {

  localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
  keys <- localKeys_1[localKeys_1$MARKER == marker, ]
  nkeys <- nrow(keys)
  if (nkeys == 0)
    return(list(results = results, processed_items = processed_items))

  # AI-043: Read the existing inference file ONCE before the per-key loop. The
  # previous design re-read the file inside every iteration AND rbind.fill'd its
  # entire content into the running `results` accumulator, so on a run that
  # iterates over (HYPO, HYPER) for the same marker the file got the
  # already-saved rows added twice — visible as N x 2 duplication in the
  # output CSV (e.g. DELTARQ_HYPO with 35292 rows instead of 17646).
  # Reading once + using the pre-loaded snapshot for the area_to_remove filter
  # keeps both behaviours correct without growing `results` across iterations.
  if (file.exists(fileNameResults) && file.info(fileNameResults)$size > 10) {
    # AI-078: polars read del CSV di resume cache, 10-30x piu' veloce di
    # utils::read.csv2 su file da 600-880 MB. Mantiene la write con
    # write.csv2 per evitare incompatibilita' di formato (polars writes
    # boolean come 'true'/'false' minuscoli, read.csv2 poi li carica come
    # character e rompe i subset logici downstream).
    # AI-061+ (2026-06-09): three polars 1.x quirks rolled into one read_csv:
    #   1. null_values="NA"          — utils::write.csv2 emits "NA" literal
    #      for missing, polars treats only "" as null on numeric columns,
    #      so without this it fails on "NA" in an f64-locked column.
    #   2. infer_schema_length large — early rows can be all zeros for
    #      INTERCEPT_PVALUE / similar, polars infers i64, then later finds a
    #      scientific-notation float (e.g. "2,52861832797769e-304") and
    #      fails to coerce to integer. Scanning more rows up-front lets it
    #      infer Float64 correctly.
    #   3. decimal_comma=TRUE        — write.csv2 uses "," as decimal sep.
    # Both quirks were exposed on ewas v32 / v33 mid-association.
    old_results_global <- unique(as.data.frame(
      polars::pl$read_csv(fileNameResults, separator = ";",
                          decimal_comma = TRUE, null_values = "NA",
                          infer_schema_length = 100000L)
    ))
    results <- plyr::rbind.fill(results, old_results_global)
  } else {
    old_results_global <- data.frame()
  }

  ssEnv_local <- get_session_info()
  tech_is_longread <- !is.null(ssEnv_local$tech) &&
                       ssEnv_local$tech %in% c("WGBS", "LONGREAD")

  for (k in seq_len(nkeys)) {
    key <- keys[k, ]
    # AI-098 (2026-06-09): symmetric tech-aware skip. Each technology has
    # exactly one canonical AREA representation; the other is no-op:
    #   - Illumina (K27/K450/K850): PROBE is canonical — literature reports
    #     probe IDs (cg00000029). POSITION would produce a duplicate
    #     coord-keyed CSV with the same numerical results → skip.
    #   - WGBS / LONGREAD: POSITION is canonical — long-reads have no
    #     "probe" concept; coordinates are the natural row identifier.
    #     PROBE pivot doesn't exist for these techs → skip.
    # This replaces the unconditional `if (AREA == "POSITION") next` which
    # was Illumina-centric and blocked all position-level inference for
    # long-reads, even when POSITION was the only meaningful unit.
    if (key$AREA == "POSITION" && !tech_is_longread) next
    if (key$AREA == "PROBE"    &&  tech_is_longread) next

    pivot_filename <- pivot_file_name_parquet(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)

    # AI-027: read via unified dispatcher. Returns NULL when neither the
    # cached parquet nor per-sample bed/bedgraph files are available,
    # which is the case run_depth_n_marker needs to skip with a warning.
    pivot_lazy <- read_pivot(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)
    if (is.null(pivot_lazy)) {
      log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
        " File not found:", pivot_filename, ".")
      association_analysis_log(cbind(prep$inference_detail, keys[k, ]),
        start_time, Sys.time(), processed_items)
      next
    }

    selected_areas_temp <- selected_areas

    # AI-061: low-memory lazy path for limma_/voom_ families. Bypass the
    # full pivot materialisation + transpose + sample_sheet merge that
    # the legacy chunked loop runs below, all of which together peak at
    # ~30× the raw matrix on a 366k-probe pivot. apply_stat_model_batch_lazy()
    # consumes the LazyFrame directly, applies the AI-043 resume filter
    # in polars, and materialises ONE R matrix only.
    is_batch_family <- grepl("^(limma|voom)_", family_test)
    if (is_batch_family) {
      area_to_remove <- character(0)
      if (nrow(old_results_global) > 0) {
        area_to_remove <- old_results_global[
          old_results_global$MARKER  == key$MARKER &
          old_results_global$FIGURE  == key$FIGURE &
          old_results_global$SUBAREA == key$SUBAREA &
          old_results_global$AREA    == key$AREA, "AREA_OF_TEST"]
      }
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Batch family '", family_test,
                "': lazy polars path (AI-061). area_to_remove=",
                length(area_to_remove))

      result_temp_local_batch <- apply_stat_model_batch_lazy(
        pivot_lazy           = pivot_lazy,
        sample_sheet         = prep$sample_names,
        family_test          = family_test,
        covariates           = prep$covariates,
        key                  = key,
        transformation_y     = prep$transformation_y,
        independent_variable = prep$independent_variable,
        area_to_remove       = area_to_remove
      )
      # AI-077: skip the save when the batch returned NULL/0 rows (resume-skip
      # case: "nothing left after resume filter"). Without this guard we
      # rewrite the entire 600-880 MB CSV identical to what's already on disk,
      # 4x per round x N round = several minutes of pure I/O waste per family.
      # The save still runs whenever new rows were appended.
      new_rows_appended <- !is.null(result_temp_local_batch) &&
                           nrow(result_temp_local_batch) > 0L
      if (new_rows_appended) {
        results <- plyr::rbind.fill(results, result_temp_local_batch)
        results <- results[, !grepl("SAMPLES_SQL_CONDITION", colnames(results)), drop = FALSE]
        association_analysis_save_results(results, fileNameResults, family_test, filter_p_value)
      }

      association_analysis_log(cbind(prep$inference_detail, keys[k, ]),
        start_time, Sys.time(), processed_items)
      if (nrow(results) != 0)
        results <- subset(results, MARKER == key$MARKER)

      # Force release of polars wrappers + apply_stat_model_batch_lazy locals
      # (y_mat ~12 GB on 367k×4k SIGNAL@PROBE + MArrayLM fit + voom weights
      # + design). R's lazy GC otherwise carries them across iterations and
      # the second batch on the same family OOMs at ~165 GB compressed
      # (limma_2 SIGNAL@PROBE, 2026-06-05). The smoke test on 2026-06-03
      # never tripped this because it ran a single batch.
      rm(result_temp_local_batch, pivot_lazy)
      gc(verbose = FALSE)
      next
    }

    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
      " Starting to read pivot:", pivot_filename, ".")
    tempDataFrame <- as.data.frame(pivot_lazy$collect())

    # AI-043: use the pre-loaded snapshot (old_results_global) for the
    # area_to_remove filter, NOT a fresh re-read of the file. Don't rbind.fill
    # old_results into the running 'results' accumulator either — that was the
    # source of cross-iteration row doubling. The file's content was already
    # folded into 'results' once, before the for-k loop opened.
    #
    # AI-062: gene names with '-' (e.g. 'HLA-A', 'ANKHD1-EIF4EBP3') are
    # rewritten to '_' by data_preparation()'s colname sanitisation, so the
    # AREA_OF_TEST that lands in the CSV uses underscores while the AREA
    # column inside the freshly-read pivot still has the dash. Without
    # normalising both sides of the %in% test, ~281 genes with dashes were
    # re-fitted on every resume run. Apply the same gsub to tempDataFrame
    # so the membership test matches the on-disk convention.
    if (nrow(old_results_global) > 0) {
      area_to_remove <- old_results_global[old_results_global$MARKER == key$MARKER &
                                            old_results_global$FIGURE == key$FIGURE &
                                            old_results_global$SUBAREA == key$SUBAREA &
                                            old_results_global$AREA == key$AREA, "AREA_OF_TEST"]
      tempDataFrame <- tempDataFrame[!(tempDataFrame$AREA %in% area_to_remove), ]
    }

    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
      " Read pivot:", pivot_filename, " with ", nrow(tempDataFrame), " rows.")
    tempDataFrame[is.na(tempDataFrame)] <- 0

    # filter by selected_areas (range or list)
    if (length(selected_areas_temp) > 0) {
      if (any(grepl(":", selected_areas_temp))) {
        selected_areas_temp <- selected_areas_temp[grepl(":", selected_areas_temp)][1]
        selected_areas_temp <- unlist(strsplit(selected_areas_temp, ":"))
        min_col <- min(as.numeric(selected_areas_temp[1]), nrow(tempDataFrame))
        max_col <- min(as.numeric(selected_areas_temp[2]), nrow(tempDataFrame))
        selected_areas_temp <- seq(from = min_col, to = max_col)
        tempDataFrame <- tempDataFrame[selected_areas_temp, ]
      } else {
        tempDataFrame <- tempDataFrame[tempDataFrame$AREA %in% selected_areas_temp, ]
      }
      if (nrow(tempDataFrame) == 0) {
        log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
          " No areas selected for the analysis! Skipped.")
      }
    }

    if (is.null(dim(tempDataFrame))) next
    if (plyr::empty(tempDataFrame) | nrow(tempDataFrame) == 0) next

    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
      " Starting to execute required test for:",
      key$MARKER, key$FIGURE, key$AREA, key$SUBAREA, ".")

    # AI-040 Fase 3: limma_<N> and voom_<N> need the WHOLE pivot at once
    # for eBayes shrinkage to be statistically meaningful. Per-chunk
    # limma estimates the prior variance from a chunk-specific subset,
    # so p-values become dependent on chunk boundaries — leaky for the
    # empirical-Bayes interpretation. Force batch families to a single
    # whole-pivot pass instead of the default chunked loop.
    chunk_size <- if (grepl("^(limma|voom)_", family_test)) {
      log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Batch family '", family_test,
                "': bypassing chunking, passing whole pivot (",
                nrow(tempDataFrame), " areas) to apply_stat_model_batch.")
      nrow(tempDataFrame)
    } else {
      ceiling(6000000 / ncol(tempDataFrame))
    }
    for (i in seq(1, nrow(tempDataFrame), by = chunk_size)) {
      chunk_indices <- i:min(i + chunk_size - 1, nrow(tempDataFrame))
      batch_df <- as.data.frame(tempDataFrame)[chunk_indices, ]
      rownames(batch_df) <- batch_df[, 1]
      batch_df <- batch_df[, -1]
      batch_df <- as.data.frame(t(batch_df))
      batch_df$Sample_ID <- rownames(batch_df)
      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
        " Transposed pivot:", pivot_filename, " with ", ncol(batch_df) - 1, " columns.")

      if (nrow(batch_df) > 1) {
        batch_df <- merge(x = prep$sample_names, y = batch_df,
                          by.x = "Sample_ID", by.y = "Sample_ID", all.x = TRUE)
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
          " Merged pivot:", pivot_filename, " with ", ncol(batch_df), " columns.")
        batch_df <- as.data.frame(batch_df)
        batch_df[is.na(batch_df)] <- 0
        batch_df <- batch_df[, -1]
        cols <- colnames(batch_df)
        batch_df <- as.data.frame(batch_df)
        if (length(colnames(batch_df)) != length(cols))
          stop("ERROR: I'm stopping here data to associate are not correct, file a bug!")
        colnames(batch_df) <- cols
        g_start <- 2 + length(prep$covariates)
        processed_items <- processed_items + ncol(batch_df) - g_start
        if (any(is.na(batch_df))) {
          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Missing values in the data frame!")
        }
        result_temp_local_batch <- apply_stat_model(
          tempDataFrame   = batch_df,
          g_start         = g_start,
          family_test     = family_test,
          covariates      = prep$covariates,
          key             = key,
          transformation_y = prep$transformation_y,
          dototal         = (length(selected_areas_temp) == 0),
          session_folder  = ssEnv$session_folder,
          prep$independent_variable,
          prep$depth_analysis,
          prep$inference_detail$samples_sql_condition,
          inference_detail = prep$inference_detail,
          ...)
        results <- plyr::rbind.fill(results, result_temp_local_batch)
        results <- results[, !grepl("SAMPLES_SQL_CONDITION", colnames(results)), drop = FALSE]
        # AI-061+ (2026-06-09): mirror the AI-077 save-guard from the
        # batch-family branch above. Only rewrite the (potentially
        # hundreds-of-MB) CSV when this chunk actually appended new
        # rows — full resume case (nothing new) should be a no-op.
        new_rows_appended_chunk <- !is.null(result_temp_local_batch) &&
                                   nrow(result_temp_local_batch) > 0L
        if (new_rows_appended_chunk) {
          association_analysis_save_results(results, fileNameResults, family_test, filter_p_value)
          n_new_rows_total <- (if (exists("n_new_rows_total")) n_new_rows_total else 0L) +
                              nrow(result_temp_local_batch)
        }
      }
    }

    association_analysis_log(cbind(prep$inference_detail, keys[k, ]),
      start_time, Sys.time(), processed_items)
    association_analysis_log(cbind(prep$inference_detail, keys[k, ]),
      start_time, Sys.time(), processed_items)
    if (nrow(results) != 0)
      results <- subset(results, MARKER == key$MARKER)
  }

  # final per-marker save (was lines 428-430)
  results <- results[, !grepl("SAMPLES_SQL_CONDITION", colnames(results)), drop = FALSE]
  association_analysis_save_results(results, fileNameResults, family_test, filter_p_value)

  list(results = results, processed_items = processed_items)
}
