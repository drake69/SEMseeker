metrics_name_collect <- function(data_frame)
{
  # NOTE: this function body is entirely commented out and currently does nothing.
  # Kept for historical reference; candidate for removal in a future cleanup.
  # # # browser
  # metrics_temp <- toupper(colnames(data_frame))
  # # read from package cache
  # if(file.exists("~/Desktop/metrics.csv"))
  #   metrics <- utils::read.csv2("~/Desktop/metrics.csv")
  # # metrics <- semseeker::metrics
  # #remove from temp existing metrics and all with PVALUE inside
  # metrics_temp <- metrics_temp[!(metrics_temp %in% semseeker::metrics_properties$Metric)]
  # metrics_temp <- metrics_temp[!(grepl("PVALUE",metrics_temp))]
  # if (length(metrics_temp)>0)
  # {
  #   metrics_temp <- data.frame("NAME"=metrics_temp)
  #   metrics_temp$IS_A_METRIC <- NA
  #   metrics_temp$SCALE_AFFECTED <- NA
  #   metrics_temp$THE_LOWER_THE_BETTER <- NA
  #   metrics <- rbind(metrics, metrics_temp)
  #   # sort by name
  #   metrics <- metrics[order(semseeker::metrics_properties$Metric),]
  #   # utils::write.csv2(metrics, file = "metrics.csv", row.names = FALSE)
  #   # usethis::use_data(metrics, overwrite = TRUE,compress = "bzip2")
  #   write.csv2(metrics, file = "~/Desktop/metrics.csv")
  #   # save to semseeker::metrics
  #   # semseeker::metrics <<- metrics
  # }
}
