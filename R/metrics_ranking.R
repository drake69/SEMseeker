metrics_ranking <- function (metric,data_frame, column_to_rank ="REBASED"){

  #
  the_lower_the_better_markers <- toupper(SEMseeker::metrics_properties[SEMseeker::metrics_properties$Higher_the_Better==FALSE,"Metric"])

  # R_SQUARED	R_SQUARED_ADJ
  the_higher_the_better_markers <-  toupper(SEMseeker::metrics_properties[SEMseeker::metrics_properties$Higher_the_Better==TRUE,"Metric"])

  # replace -Inf with 1E-ssEnv$plot_resolution
  data_frame[data_frame == -Inf] <- 1E-300
  # replace Inf with 1E300
  data_frame[data_frame == Inf] <- 1E300

  if(grepl("PVALUE",metric))
    the_lower_the_better_markers <- c(the_lower_the_better_markers, metric)

  # metric <- toupper(metric)
  # # data_frame[,column_to_rank] <- abs(data_frame[,column_to_rank])
  # if(metric %in% the_lower_the_better_markers){
  #   data_frame <- data_frame[order(data_frame[,column_to_rank], decreasing = T),]
  # }else if(metric %in% the_higher_the_better_markers){
  #   data_frame <- data_frame[order(data_frame[,column_to_rank], decreasing = F),]
  # }else{
  #   log_event("ERROR: Metric not found")
  #   stop()
  # }
  #
  #
  # values_to_rank <- unique(data_frame[,column_to_rank])
  # values_to_rank <- data.frame("VALUE"=values_to_rank)
  # # create a new column with the rank
  # # values_to_rank$SCORE <- 1:nrow(values_to_rank)
  #
  # values_to_rank$SCORE <- cut(values_to_rank$VALUE, breaks = 100, right = FALSE, labels = FALSE)
  # # merge the rank with the data frame
  # data_frame <- merge(data_frame, values_to_rank, by.x=column_to_rank, by.y="VALUE", all.x=TRUE)

  data_frame$METRIC <- metric

  if(max(data_frame[,column_to_rank]) == min(data_frame[,column_to_rank])){
    data_frame$SCORE <- 1
    return(data_frame)
  }


  if(metric %in% the_lower_the_better_markers){
    data_frame$SCORE <- normalize_minimize(data_frame[,column_to_rank])
  }else if(metric %in% the_higher_the_better_markers){
    data_frame$SCORE <- normalize_maximize(data_frame[,column_to_rank])
  }else{
    log_event("ERROR: Metric not found")
    stop()
  }

  if(max(data_frame$SCORE)>1)
  {

  }
  data_frame$SCORE <- rank(data_frame$SCORE, ties.method = "min")
  # # add the marker and figure to the data frame
  # scores <- rbind(scores, )

  return(data_frame)
}


# Normalisation helpers used by metrics_ranking():
#   normalize_minimize: maps x to [0,1] where the MAXIMUM of x → 0 and
#     the MINIMUM → 1 (lower-is-better metrics get a higher normalised score).
#   normalize_maximize: maps x to [0,1] where the MINIMUM of x → 0 and
#     the MAXIMUM → 1 (higher-is-better metrics get a higher normalised score).
normalize_minimize <- function(x) {(max(x) - x) / (max(x) - min(x))}
normalize_maximize <- function(x) {(x - min(x)) / (max(x) - min(x))}
