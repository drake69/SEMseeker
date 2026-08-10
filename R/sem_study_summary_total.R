sem_study_summary_total <- function()
{

  ssEnv <- core_get_session_info()
  study_summary <- sem_study_summary_get()

  # ssEnv <- core_get_session_info()
  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- subset(keys, AREA=="POSITION")
  if(nrow(keys)==0)
    return()

  # AI-083: temp_result MUST be a local binding initialised here. The previous
  # code used exists("temp_result"), which resolves through the enclosing
  # environments up to globalenv(): a leftover object of that name in an
  # interactive session (or restored from a .RData) made the first iteration
  # take the merge branch against an unrelated data frame, so every burden
  # column landed as NA while PROBES_COUNT (assigned unconditionally below)
  # stayed populated — exactly the reported symptom.
  temp_result <- NULL

  for ( k in seq_len(nrow(keys)))
  {
    # k <- 1
    key <- keys[k,]
    marker <- key$MARKER
    figure <- key$FIGURE
    area <- key$AREA
    subarea <- key$SUBAREA
    # AI-027: read via unified dispatcher. NULL means no pivot AND no
    # streaming-merge source — skip the key with the same semantics as
    # the previous file.exists() guard.
    pivot <- io_read_pivot(marker, figure, area, subarea)
    if (is.null(pivot))
      next
    # remove CHR START END columns
    pivot <- pivot$drop(c("CHR", "START", "END"))
    # sum per columns
    if(key$DISCRETE)
      pivot <- pivot$sum()$with_columns(polars::pl$col("*"))
    else
      pivot <- pivot$mean()$with_columns(polars::pl$col("*"))

    pivot <- as.data.frame(t(as.data.frame(pivot$collect())))
    combined_key <- paste0(marker,"_",figure)
    colnames(pivot) <- combined_key
    pivot$Sample_ID <- row.names(pivot)

    if(is.null(temp_result))
      temp_result <- pivot
    else
    {
      temp_result <- temp_result[, !(colnames(temp_result) == combined_key), drop = FALSE]
      temp_result <- merge(temp_result, pivot, by="Sample_ID", all=TRUE)
    }
  }
  if (is.null(temp_result))
    return()
  temp_result[is.na(temp_result)] <- 0
  # remove from summary all column excet Sample_ID from temp_result
  col_temp <- colnames(temp_result)[!(colnames(temp_result) == "Sample_ID")]
  study_summary <- study_summary[, !(colnames(study_summary) %in% col_temp)]

  # AI-083: the burden merge is name-based. If the two sides share no
  # Sample_ID the all.x=TRUE merge silently yields a sample sheet whose burden
  # columns are 100% NA (PROBES_COUNT below still gets filled), which then
  # crashes every depth=1 inference downstream with "data are not the same
  # size". Fail here, naming both sides, instead of writing that file.
  ids_summary <- as.character(study_summary$Sample_ID)
  ids_burden  <- as.character(temp_result$Sample_ID)
  matched     <- intersect(ids_summary, ids_burden)
  if (length(matched) == 0) {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Burden aggregation matches no sample: sample sheet ids [",
              paste(utils::head(ids_summary, 4), collapse = ", "),
              "] vs pivot columns [",
              paste(utils::head(ids_burden, 4), collapse = ", "), "].")
    stop("sem_study_summary_total(): no Sample_ID in common between the sample ",
         "sheet and the POSITION pivots — the burden columns would all be NA. ",
         "Sample sheet ids: ", paste(utils::head(ids_summary, 4), collapse = ", "),
         " | pivot columns: ", paste(utils::head(ids_burden, 4), collapse = ", "), ".")
  }
  missing_ids <- setdiff(unique(ids_summary), ids_burden)
  if (length(missing_ids) > 0)
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", length(missing_ids), " sample(s) have no burden in the POSITION pivots: ",
              paste(utils::head(missing_ids, 10), collapse = ", "),
              if (length(missing_ids) > 10) ", ..." else "")

  study_summary <- merge(study_summary, temp_result, by="Sample_ID", all.x=TRUE)
  study_summary$PROBES_COUNT <- ssEnv$probes_count
  # io_file_path_build() uppercases via core_name_cleaning() — on-disk name is
  # SAMPLE_SHEET_RESULT.csv. Linux ext4 is case-sensitive; readers that
  # hard-code the path must use the uppercase form. See io_file_path_build().
  summary_file <- io_file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  utils::write.csv2(study_summary,summary_file, row.names = FALSE)
}
