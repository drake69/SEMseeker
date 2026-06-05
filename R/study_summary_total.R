study_summary_total <- function()
{

  ssEnv <- get_session_info()
  study_summary <- study_summary_get()

  # ssEnv <- get_session_info()
  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- subset(keys, AREA=="POSITION")
  if(nrow(keys)==0)
    return()
  for ( k in 1:nrow(keys))
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
    pivot <- read_pivot(marker, figure, area, subarea)
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

    if(!exists("temp_result"))
      temp_result <- pivot
    else
    {
      temp_result <- temp_result[, !(colnames(temp_result) == combined_key)]
      temp_result <- merge(temp_result, pivot, by="Sample_ID", all=TRUE)
    }
  }
  if (!exists("temp_result"))
    return()
  temp_result[is.na(temp_result)] <- 0
  # remove from summary all column excet Sample_ID from temp_result
  col_temp <- colnames(temp_result)[!(colnames(temp_result) == "Sample_ID")]
  study_summary <- study_summary[, !(colnames(study_summary) %in% col_temp)]
  study_summary <- merge(study_summary, temp_result, by="Sample_ID", all.x=TRUE)
  study_summary$PROBES_COUNT <- ssEnv$probes_count
  summary_file <- file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  utils::write.csv2(study_summary,summary_file, row.names = FALSE)
}
