assoc_analysis_log <- function(inference_detail, start_time, end_time, processed_items)
{
  ssEnv <- core_get_session_info()
  log_folder <- ssEnv$session_folder
  association_file <- paste0(log_folder, "/association_analysis.csv")
  inference_detail$node_name <- as.character(Sys.info()["nodename"])
  inference_detail$session_id <- ssEnv$session_id
  # inference_detail$start_time <- as.character(format(start_time, "%a %b %d %X %Y"))
  # inference_detail$end_time <- as.character(format(end_time, "%a %b %d %X %Y"))
  inference_detail$processed_time <- as.numeric(difftime(end_time, start_time, units  =  "mins"))
  inference_detail$processed_items <- processed_items
  if(!file.exists(association_file))
  {
    utils::write.csv2(inference_detail, file  =  association_file, row.names  =  FALSE)
  } else
  {
    # convert all columns of inference detail as character
    tryCatch({
      association_file_data <- utils::read.csv2(association_file, header  =  TRUE, stringsAsFactors  =  FALSE)
      association_file_data <- plyr::rbind.fill(inference_detail, association_file_data)
      utils::write.csv2(association_file_data, file  =  association_file, row.names  =  FALSE)
    }, error = function(e) {
      # print("ERROR: I'm stopping here, data to associate are not correct, file a bug!")
    })
  }

}
