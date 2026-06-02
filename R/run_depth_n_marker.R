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

  for (k in seq_len(nkeys)) {
    key <- keys[k, ]
    if (key$AREA == "POSITION") next
    pivot_filename <- pivot_file_name_parquet(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)

    if (!file.exists(pivot_filename)) {
      log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
        " File not found:", pivot_filename, ".")
      association_analysis_log(cbind(prep$inference_detail, keys[k, ]),
        start_time, Sys.time(), processed_items)
      next
    }

    selected_areas_temp <- selected_areas
    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),
      " Starting to read pivot:", pivot_filename, ".")
    tempDataFrame <- as.data.frame(polars::pl$read_parquet(pivot_filename))

    if (file.exists(fileNameResults) && file.info(fileNameResults)$size > 10) {
      old_results <- unique(utils::read.csv2(fileNameResults, header = TRUE))
      area_to_remove <- old_results[old_results$MARKER == key$MARKER &
                                     old_results$FIGURE == key$FIGURE &
                                     old_results$SUBAREA == key$SUBAREA &
                                     old_results$AREA == key$AREA, "AREA_OF_TEST"]
      tempDataFrame <- tempDataFrame[!(tempDataFrame$AREA %in% area_to_remove), ]
      results <- plyr::rbind.fill(results, old_results)
      rm(old_results)
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

    chunk_size <- ceiling(6000000 / ncol(tempDataFrame))
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
          ...)
        results <- plyr::rbind.fill(results, result_temp_local_batch)
        results <- results[, !grepl("SAMPLES_SQL_CONDITION", colnames(results)), drop = FALSE]
      }
      association_analysis_save_results(results, fileNameResults, family_test, filter_p_value)
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
