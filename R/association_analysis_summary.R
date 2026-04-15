association_analysis_summary <- function(inference_details,destination_folder="", result_folder="", ...)
{
  if(result_folder!="")
    ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  else
    ssEnv <- get_session_info()
  association_data <- association_data_extractor(inference_details, destination_folder, result_folder, ...)


  available_metrics <- toupper(SEMseeker::metrics_properties[,"Metric"])

  if(any(grepl("PVALUE",colnames(association_data))))
    available_metrics <- c(available_metrics, colnames(association_data)[grepl("PVALUE",colnames(association_data))])

  # remove not numeric columns metrics
  metrics_to_remove <- colnames(association_data[,!sapply(association_data, is.numeric)])
  available_metrics <- available_metrics[!(available_metrics %in% metrics_to_remove)]

  available_metrics <- available_metrics[available_metrics %in% colnames(association_data)]

  #  create a summary table for the association analysis grouping by AREA,SUBAREA,MARKER,FIGURE and SAMPLES_SQL_CONDITION if exists
  if(any("SAMPLES_SQL_CONDITION" %in% colnames(association_data))) {
    summary_table <- association_data %>%
      dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$SAMPLES_SQL_CONDITION) %>%
      dplyr::summarise(dplyr::across(available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  } else {
    summary_table <- association_data %>%
      dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE) %>%
      dplyr::summarise(dplyr::across(available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  }

  summary_table <- as.data.frame(summary_table)
  return(summary_table)
}
