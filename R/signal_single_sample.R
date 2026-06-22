#' signal_single_sample
#'
#' @param values signal vaòues
#' @param sample_detail detais of sample
#' @param probe_features annotation probe
#'
#' @return signal mean
#'
signal_single_sample <- function(values,sample_detail,probe_features)
{
  ssEnv <- get_session_info()

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData, c(sample_detail$Sample_Group ,paste0("SIGNAL","_", "MEAN", sep = "")))
  # DEBUG (2026-06-09): right before the data.frame() that has
  # been exploding with "arguments imply differing number of rows" since
  # v35. Inspect:
  #   length(values), nrow(probe_features), class(values)
  #   head(names(values)), head(probe_features$PROBE)
  #   is.null(names(values)), all(probe_features$PROBE %in% names(values))
 
  signal_values_annotated <- data.frame(as.data.frame(probe_features), "VALUE" = values, row.names = probe_features$PROBE)[, c("CHR", "START", "END","VALUE")]
  dump_sample_as_bed_file(
    data_to_dump = signal_values_annotated,
    fileName = file_path_build(baseFolder =  folder_to_save, detailsFilename =  c(sample_detail$Sample_ID,"SIGNAL","MEAN"), extension = "bedgraph", add_gz=TRUE)
  )
}
