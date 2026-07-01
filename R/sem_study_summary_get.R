sem_study_summary_get <- function(sql_sample_selection="")
{
  ssEnv <- core_get_session_info()
  # io_file_path_build() uppercases via core_name_cleaning() — the on-disk file is
  # SAMPLE_SHEET_RESULT.csv, written by sem_study_summary_total(). Going through
  # io_file_path_build() here guarantees the same path is resolved on every OS
  # regardless of file-system case sensitivity.
  summary_file <- io_file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv")
  if(!file.exists(summary_file))
    summary_file <- io_file_path_build( ssEnv$result_folderData, "1_SAMPLE_SHEET_ORIGINAL","csv")

  if(!file.exists(summary_file))
  {
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No sample sheet found in ", ssEnv$result_folderData)
    return()
  }

  study_summary <-   utils::read.csv2(summary_file)
  study_summary <- assoc_filter_sql(sql_sample_selection, study_summary)

  reference <- study_summary[study_summary$Sample_Group == "Reference",]
  no_reference <- study_summary[study_summary$Sample_Group != "Reference",]
  reference <- reference[!(reference$Sample_ID %in% no_reference$Sample_ID), ]

  study_summary <- rbind(no_reference, reference)

  if(nrow(reference)==0)
    study_summary <- no_reference

  return(study_summary)
}
