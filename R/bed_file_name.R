bed_file_name <- function(sample_id,sample_group,marker,figure)
{
  ssEnv <- get_session_info()
  if (is.na(sample_group) || sample_group=="")
  {
    stop("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " bed_file_name: sample_group is empty")
  }

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_group),paste0(marker,"_",figure)))
  bed_ext <- unique(ssEnv$keys_markers_figures_default[ ssEnv$keys_markers_figures_default$MARKER==marker & ssEnv$keys_markers_figures_default$FIGURE==figure, "EXT"])
  if(length(bed_ext)==0){
    stop("ERROR: bed_file_name: bed_ext is empty")
  }
  file_name <- file_path_build(folder_to_save,c(sample_id,marker,figure),bed_ext, add_gz=TRUE)
  return(file_name)
}
