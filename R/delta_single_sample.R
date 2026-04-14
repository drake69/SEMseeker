#' delta_single_sample
#'
#' @param values data.frame of methylation values with columns CHR, START, END and signal value in column 4
#' @param thresholds data.frame of signal thresholds (from signal_range_values) with columns CHR, START, END,
#'   signal_superior_thresholds, signal_inferior_thresholds, signal_median_values
#' @param sample_detail named list/row with at least Sample_ID and Sample_Group fields
#' @return invisibly NULL; results are written as bedgraph.gz files under the session data folder
#'
delta_single_sample <- function(values, thresholds, sample_detail) {

  ssEnv <- get_session_info()

  # Polars inner join on (CHR, START, END) — replaces sort-then-positional-zip.
  # join_values_to_thresholds() handles type coercion and is shared with
  # mutations_get() and deltar_single_sample().
  joined <- join_values_to_thresholds(values, thresholds)

  if (nrow(joined) == 0L) {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " [delta_single_sample] No overlapping positions — skipping delta",
      " for sample=", sample_detail$Sample_ID)
    return(invisible(NULL))
  }

  high_thresholds <- joined$signal_superior_thresholds
  low_thresholds  <- joined$signal_inferior_thresholds

  if (any(high_thresholds < low_thresholds))
    stop("ERROR: I'm stopping here — some high thresholds have values less than low thresholds!")

  ### get deltas HYPER ###########################################################
  deltas_hyper <- data.frame(
    CHR   = joined$CHR,
    START = joined$START,
    END   = joined$END,
    DELTA = joined$VALUE - high_thresholds,
    stringsAsFactors = FALSE
  )
  deltas_hyper <- sort_by_chr_and_start(deltas_hyper)
  deltas_hyper <- subset(deltas_hyper, deltas_hyper$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,
    c(as.character(sample_detail$Sample_Group), "DELTAS_HYPER"))
  dump_sample_as_bed_file(
    data_to_dump = deltas_hyper,
    fileName     = file_path_build(folder_to_save,
                     c(as.character(sample_detail$Sample_ID), "DELTAS", "HYPER"),
                     "bedgraph", add_gz = TRUE)
  )

  ### get deltas HYPO ############################################################
  deltas_hypo <- data.frame(
    CHR   = joined$CHR,
    START = joined$START,
    END   = joined$END,
    DELTA = low_thresholds - joined$VALUE,
    stringsAsFactors = FALSE
  )
  deltas_hypo <- sort_by_chr_and_start(deltas_hypo)
  deltas_hypo <- subset(deltas_hypo, deltas_hypo$DELTA > 0)[, c("CHR", "START", "END", "DELTA")]

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,
    c(as.character(sample_detail$Sample_Group), "DELTAS_HYPO"))
  dump_sample_as_bed_file(
    data_to_dump = deltas_hypo,
    fileName     = file_path_build(folder_to_save,
                     c(as.character(sample_detail$Sample_ID), "DELTAS", "HYPO"),
                     "bedgraph", add_gz = TRUE)
  )

  deltas_to_check <- c(deltas_hypo$DELTA, deltas_hyper$DELTA)
  if (length(deltas_to_check) > 0)
    if (min(deltas_to_check) < 0)
      stop("ERROR: I'm stopping here — deltas have negative values!")
}
