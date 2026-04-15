pivot_file_name <- function(marker,figure,area,subarea,add_gz=TRUE)
{
  ssEnv <- get_session_info()
  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  pivot_subfolder <- dir_check_and_create(reportFolder, marker)
  # C-06: append genome_build suffix (e.g. "_hg19") for belt-and-suspenders provenance
  genome_suffix <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    ssEnv$genome_build else "hg19"
  pivot_base <- paste0(marker,"_", figure,"_",area,"_",subarea,"_",genome_suffix, sep = "")
  pivot_file_name <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_base,extension =  ".csv" ,add_gz=add_gz)
  return(pivot_file_name)
}

pivot_file_name_parquet <- function(marker,figure,area,subarea)
{
  ssEnv <- get_session_info()
  reportFolder <- dir_check_and_create(ssEnv$result_folderData,"Pivots")
  pivot_subfolder <- dir_check_and_create(reportFolder, marker)
  # C-06: append genome_build suffix (e.g. "_hg19") for belt-and-suspenders provenance
  genome_suffix <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    ssEnv$genome_build else "hg19"
  pivot_base <- paste0(marker,"_", figure,"_",area,"_",subarea,"_",genome_suffix, sep = "")
  pivot_file_name <- file_path_build(baseFolder =  pivot_subfolder,detailsFilename =  pivot_base,extension =  ".parquet" ,add_gz=FALSE)
  return(pivot_file_name)
}
