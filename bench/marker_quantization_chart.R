#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
# NOTE: currently internal — could be re-exported once @examples are added
marker_quantization_chart <- function(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{
  # result_folder, maxResources = 90, parallel_strategy  = "multisession", ...
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)
  # ssEnv <- get_session_info()

  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- keys[keys$AREA == "POSITION",]
  if(nrow(keys) == 0)
  {
    log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"), " For this analisys the POSITION area is required")
    close_env()
    return()
  }

  keys <- keys[keys$MARKER != "DELTAS",]
  keys <- keys[keys$MARKER != "DELTAR",]
  keys <- keys[keys$MARKER != "LESIONS",]
  keys <- keys[keys$MARKER != "SIGNAL",]

  keys <- keys[complete.cases(keys),]
  nkeys <- nrow(keys)

  #
  #

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:nkeys)
  else
    progress_bar <- ""


  to_export <- c("keys","result_folderPivot","ssEnv","file_path_build","progress_bar", "progression_index","progression","progressor_uuid","owner_session_uuid","trace")
  result_temp <- data.frame()
  scores <- data.frame()


  # result_temp <- foreach::foreach(k = 1:nkeys, .combine =  plyr::rbind.fill, .export = to_export) %dorng%
  for (k in 1:nkeys)
  {

    # k <- 1
    key <- keys [k,]
    log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), "Processing key: ", key$MARKER," ", key$FIGURE," ", key$AREA," ", key$SUBAREA)
    if(is.na(key$MARKER))
      next

    if(key$MARKER=="DELTARP" | key$MARKER=="DELTARQ")
      original_marker <- "DELTAR"
    else
      original_marker <- "DELTAS"


    fname <- pivot_file_name_parquet(original_marker, key$FIGURE, key$AREA,key$SUBAREA)
    if(!file.exists(fname))
    {
      log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," File not found: ", fname)
      next
    }
    original <-  polars::pl$read_parquet(fname)

    fname <- pivot_file_name_parquet(key$MARKER, key$FIGURE, key$AREA,key$SUBAREA)
    if(!file.exists(fname))
    {
      log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y")," File not found: ", fname)
      next
    }
    quantized <- polars::pl$read_parquet(fname)

    # remove zeros from original and quantized
    # remove first 3 columns and transform as vector
    quantized <- as.vector(as.matrix(as.data.frame(quantized$drop("CHR","START","END"))))
    # remove na
    quantized <- quantized[!is.na(quantized)]

    original <- as.vector(as.matrix(as.data.frame(original$drop("CHR","START","END"))))
    original <- original[!is.na(original)]

    if (key$MARKER=="MUTATIONS")
      quantized <- c(quantized,0,2)
    original_range <- range(original)
    quantized_range <- range(quantized)
    original_to_plot <- (original * (quantized_range[2] - quantized_range[1]) / (original_range[2] - original_range[1])) + (quantized_range[1] - original_range[1])
    quantized_to_plot <- (quantized * (original_range[2] - original_range[1]) / (quantized_range[2] - quantized_range[1])) + (original_range[1] - quantized_range[1])
    # Define number of intervals
    num_intervals_original <- quantized_range[2]  # Change this to desiblue number of intervals

    chart_folder <- dir_check_and_create(ssEnv$result_folderChart, "MARKERS_DISTRIBUTION")

    filename <- file_path_build( chart_folder ,c(quantized_range[2],original_marker,key$MARKER, key$FIGURE, key$AREA,key$SUBAREA, "dens"),"png")
    grDevices::png(filename, width = 9, height = 9, units="in", res = as.numeric(ssEnv$plot_resolution_ppi))
    plot(density(original), col = "skyblue", main = paste0(original_marker, " Vs. ", key$MARKER ," Density Plot"))
    lines(density(quantized_to_plot), col = "blue")
    legend("bottom", legend = c(original_marker, key$MARKER),  col = c("skyblue", "blue"), lty = 1, cex = 0.8, bty = "n")
    dev.off()

    # plot two istograms with log10 y axis scale
    filename <- file_path_build( chart_folder ,c(quantized_range[2],original_marker,key$MARKER, key$FIGURE, key$AREA,key$SUBAREA, "hist"),"png")
    grDevices::png(filename, width = 9, height = 9, units="in", res = as.numeric(ssEnv$plot_resolution_ppi))
    par(mfrow=c(2,1))
    hist(original, col = "skyblue", xlab = "Frequency", main = paste0( original_marker," Histogram"))
    hist(quantized_to_plot, col = "blue", xlab = "Frequency", main = paste0(key$MARKER, " Histogram"))
    par(mfrow=c(1,1))
    dev.off()


    if(ssEnv$showprogress)
      progress_bar(sprintf("I'm doing plot for the quantization process."))
    return(1)
  }

}



