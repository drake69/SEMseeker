#' analyze_single_sample
#'
#' @param values values of methylation
#' @param thresholds threshold to use for comparison
#' @param sample_detail details of the sample to analyze
#' @param figure which figure of sample will be analysed: HYPO or HYPER
#' @return list of lesion count and probe_features count
#'
#' @importFrom doRNG %dorng%
analyze_single_sample <- function(values, thresholds, figure, sample_detail) {

  ssEnv <- get_session_info()
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  "analyze_single_sample:", ssEnv$result_folderData)

  mutation_annotated_sorted <- mutations_get(values, figure,thresholds, sample_detail$Sample_ID)
  mutation_annotated_sorted_to_save <- subset(mutation_annotated_sorted, mutation_annotated_sorted$MUTATIONS == 1)[, c("CHR", "START", "END")]
  if (nrow(mutation_annotated_sorted_to_save) != 0)
    mutation_annotated_sorted_to_save$VALUE <- 1

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("MUTATIONS","_", figure, sep = "")))
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  "analyze_single_sample:",folder_to_save)
  dump_sample_as_bed_file(
    data_to_dump = mutation_annotated_sorted_to_save,
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"MUTATIONS",figure),"bed", add_gz=TRUE)
  )

  lesionWeighted <- lesions_get(grouping_column = "CHR", mutation_annotated_sorted = mutation_annotated_sorted)
  if(nrow(lesionWeighted) != 0)
    lesionWeighted$VALUE <- 1

  folder_to_save <- dir_check_and_create(ssEnv$result_folderData,c(as.character(sample_detail$Sample_Group),paste0("LESIONS","_", figure, sep = "")))
  dump_sample_as_bed_file(
    data_to_dump = lesionWeighted,
    fileName = file_path_build(folder_to_save,c(sample_detail$Sample_ID,"LESIONS",figure),"bed", add_gz=TRUE)
  )

}
