available_metrics <- function(inference_details,result_folder, ...)
{
  
  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  localKeys <- ssEnv$keys_markers_figures
  markers <- unique(localKeys$MARKER)
  for(z in 1:nrow(inference_details))
  {
    start_time <- Sys.time()
    inference_detail <- inference_details[z,]
    for (a in seq_along(markers) )
    {
      if (exists("results_inference"))
        rm(results_inference)

      keys <- localKeys[localKeys$MARKER==markers[a],]
      keys <- unique(keys)
      fileNameResults <- inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference)
      if (file.exists(fileNameResults))
      {
        results_inference <- utils::read.csv2(fileNameResults, header  =  TRUE, nrows = 1)
      }
      colnames_ass <- toupper(colnames(results_inference))
      pvalue_col <- colnames_ass[grepl("PVALUE", colnames_ass)]
      metrics_col <- colnames_ass[ colnames_ass %in% ssEnv$model_metrics]

      metrics_col <- c(metrics_col, pvalue_col)
      # sort
      metrics_col <- sort(metrics_col)
      # replace . and space with underscore
      metrics_col <- gsub("[.]","_",metrics_col)
      metrics_col <- gsub(" ","_",metrics_col)
      return(metrics_col)
    }
  }
}
