#--- signal_range_values ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#' calculate the range of signal values to define the outlier
#' @param populationMatrix matrix of methylation for the population under calculation (probes × samples)
#' @param batch_id character string identifying the batch; used to name the cached parquet output file
#' @param probe_features data.frame of probe annotations (from probe_features_get) with columns PROBE, CHR, START, END
#'
#' @return data.frame of per-probe thresholds with columns signal_inferior_thresholds,
#'   signal_superior_thresholds, signal_median_values, iqr, q1, q3, PROBE, CHR, START, END
#' @importFrom doRNG %dorng%

signal_range_values <- function(populationMatrix, batch_id, probe_features) {


  ssEnv <- get_session_info()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting signal thresholds calculation.")
  thresholds_file_name <- file_path_build(ssEnv$result_folderData ,c(batch_id, "signal_thresholds"),"parquet")
  if(file.exists(thresholds_file_name))
  {
    result <- as.data.frame(polars::pl$read_parquet(thresholds_file_name))
    return(result)
  }
  if (sum(is.na(populationMatrix)) > 0)
  {
    msg <- paste0("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), " There are missing values in the population matrix, apply the parameter inpute or remove the missing values.")
    log_event(msg)
    stop(msg)
  }

  # populationMatrix <- result
  # populationMatrixDim <- dim(populationMatrix)
  populationMatrix <- as.data.frame(populationMatrix)
  populationMatrix <- populationMatrix[, !(colnames(populationMatrix) %in% "PROBE")]
  # min_values <- apply(populationMatrix, 1, min, na.rm=TRUE)
  # max_values <- apply(populationMatrix, 1, max, na.rm=TRUE)
  # q1 <- apply(populationMatrix, 1, quantile, probs = 0.25, na.rm=TRUE)
  # q3 <- apply(populationMatrix, 1, quantile, probs = 0.75, na.rm=TRUE)
  # median_values <- apply(populationMatrix, 1, median, na.rm=TRUE)
  # iqr_values <- q3 - q1
  # rm(populationMatrix)
  # signal_superior_thresholds <- q3 + 3 * iqr_values
  # signal_inferior_thresholds[signal_superior_thresholds>max_values] <- max_values
  # signal_inferior_thresholds <- q1 - 3 * iqr_values
  # signal_inferior_thresholds[signal_inferior_thresholds<min_values] <- min_values

  row.names(populationMatrix) <- rownames(populationMatrix)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting signal thresholds calculation.")
  export = c("progress_bar","progression_index", "progression", "progressor_uuid", "owner_session_uuid", "trace","populationMatrix","ssEnv")
  r <- 1

  # if(ncol(populationMatrix) < 1000)
  # {
  #   # for(r in 1:1000)
  #   result <- foreach::foreach(r = 1:nrow(populationMatrix), .combine = "rbind", .export = export) %dorng%
  #     {
  #       signal_row <- populationMatrix[r,]
  #       signal_row <- as.vector(t(signal_row))
  #       temp <- stats::quantile(signal_row, c(0.25,0.5,0.75), na.rm=TRUE)
  #       signalQ1Values <-  temp[1]
  #       signalQ3Values <- temp[3]
  #       signal_median_values <- temp[2]
  #       signalValuesIQR <- stats::IQR(signal_row)
  #
  #       signal_inferior_thresholds <- max((signalQ1Values - (ssEnv$iqrTimes * signalValuesIQR)), min(signal_row, na.rm = TRUE), na.rm=TRUE)
  #       signal_superior_thresholds <- min((signalQ3Values + (ssEnv$iqrTimes * signalValuesIQR)), max(signal_row, na.rm = TRUE), na.rm=TRUE)
  #
  #       temp_result <- data.frame(
  #         "signal_inferior_thresholds"= signal_inferior_thresholds,
  #         "signal_superior_thresholds"= signal_superior_thresholds,
  #         "signal_median_values"= signal_median_values,
  #         "iqr" = signalValuesIQR,
  #         "q1"= signalQ1Values,
  #         "q3"= signalQ3Values)
  #       row.names(temp_result) <- row.names(signal_row)
  #       # temp_result$PROBE <- row.names(signal_row)
  #       # colnames(temp_result) <- c("signal_inferior_thresholds","signal_superior_thresholds","signal_median_values")
  #       if(ssEnv$showprogress)
  #         progress_bar(sprintf("%s",names(signal_row)))
  #       temp_result
  #     }
  #
  #
  # }
  # else


  {
    chunk_size <- 10000  # Define a chunk size
    result <- data.frame()
    if(ssEnv$showprogress)
      progress_bar <- progressr::progressor(along = seq(1, nrow(populationMatrix), by = chunk_size))
    for (i in seq(1, nrow(populationMatrix), by = chunk_size)) {

      chunk_indices <- i:min(i + chunk_size - 1, nrow(populationMatrix))
      th <- future.apply::future_apply(populationMatrix[chunk_indices, ], 1 ,  get_th <- function(signal_row)
      {
        signal_row <- as.vector((signal_row))
        temp <- stats::quantile(signal_row, c(0.25,0.5,0.75), na.rm = TRUE)
        signalQ1Values <-  temp[1]
        signal_median_values <- temp[2]
        signalQ3Values <- temp[3]
        signalValuesIQR <- stats::IQR(signal_row)

        signal_inferior_thresholds <- max((signalQ1Values - (as.numeric(ssEnv$iqrTimes) * signalValuesIQR)), min(signal_row, na.rm = TRUE), na.rm = TRUE)
        signal_superior_thresholds <- min((signalQ3Values + (as.numeric(ssEnv$iqrTimes) * signalValuesIQR)), max(signal_row, na.rm = TRUE), na.rm = TRUE)

        temp_result <- c(
          "signal_inferior_thresholds"= signal_inferior_thresholds,
          "signal_superior_thresholds"= signal_superior_thresholds,
          "signal_median_values"= signal_median_values,
          "iqr" = signalValuesIQR,
          "q1" = signalQ1Values,
          "q3" = signalQ3Values)
        names(temp_result) <- names(signal_row)
        # temp_result$PROBE <- names(signal_row)
        temp_result
      }, future.chunk.size = 1000)
      #
      result <- rbind(result, as.data.frame(t(th)))
      rm(th)
      if(ssEnv$showprogress)
        progress_bar(sprintf("Done %s positions",i))
    }
    gc()
  }

  #
  colnames(result) <- c("signal_inferior_thresholds","signal_superior_thresholds", "signal_median_values","iqr","q1","q3")
  result$PROBE <- rownames(populationMatrix)
  if(nrow(result) != nrow(populationMatrix))
    stop("I'M STOPPING HERE, No thresholds defined for the population.")

  result <- polars::as_polars_df(result)$lazy()
  probe_features <- polars::as_polars_df(probe_features)$lazy()
  # rename AREA as PROBE
  result <- probe_features$join(
    result,
    on = c("PROBE"),
    how = "inner"
  )
  result <- result$drop(c("PROBE_WHOLE"))
  result <- result$sort(c("CHR", "START","END"), descending = FALSE)
  result <- result$collect()
  result$write_parquet(thresholds_file_name)

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Thresholds defined for: ", nrow(result), " probe_features.")
  gc()
  return(as.data.frame(result))
}
