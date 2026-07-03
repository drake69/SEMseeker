assoc_analysis_summary <- function(inference_details,destination_folder="", result_folder="", ...)
{
  if(result_folder!="")
    ssEnv <- core_init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  else
    ssEnv <- core_get_session_info()
  association_data <- assoc_data_extractor(inference_details, destination_folder, result_folder, ...)


  sem_available_metrics <- toupper(SEMseeker::metrics_properties[,"Metric"])

  if(any(grepl("PVALUE",colnames(association_data))))
    sem_available_metrics <- c(sem_available_metrics, colnames(association_data)[grepl("PVALUE",colnames(association_data))])

  # remove not numeric columns metrics
  metrics_to_remove <- colnames(association_data[,!vapply(association_data, is.numeric, logical(1))])
  sem_available_metrics <- sem_available_metrics[!(sem_available_metrics %in% metrics_to_remove)]

  sem_available_metrics <- sem_available_metrics[sem_available_metrics %in% colnames(association_data)]

  #  create a summary table for the association analysis grouping by AREA,SUBAREA,MARKER,FIGURE and SAMPLES_SQL_CONDITION if exists
  if(any("SAMPLES_SQL_CONDITION" %in% colnames(association_data))) {
    summary_table <- association_data %>%
      dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$SAMPLES_SQL_CONDITION) %>%
      dplyr::summarise(dplyr::across(sem_available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  } else {
    summary_table <- association_data %>%
      dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE) %>%
      dplyr::summarise(dplyr::across(sem_available_metrics, list(
        max= ~max(., na.rm=TRUE),
        min= ~min(., na.rm=TRUE),
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        count_below_0.05 = ~sum(. < 0.05, na.rm = TRUE))))
  }

  summary_table <- as.data.frame(summary_table)
  return(summary_table)
}
