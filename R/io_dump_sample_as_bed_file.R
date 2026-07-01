#' given data and colnames dump as bed file
#'
#' @param data_to_dump data frame to dump into bed file with CHR, START, END
#' @param fileName name of the file to save data in
#'
#' @return nothing
#'

io_dump_sample_as_bed_file <- function(data_to_dump, fileName) {

  ssEnv <- core_get_session_info()
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " io_dump_sample_as_bed_file ssEnv:", length(ssEnv))
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " io_dump_sample_as_bed_file:", ssEnv$result_folderData)

  if (!plyr::empty(data_to_dump)) {
    data_to_dump[, "CHR"] <- anno_normalize_chr(data_to_dump[, "CHR"], "output")
  }

  # bed coordinate must start from zero!
  data_to_dump$START <- as.numeric(data_to_dump$START)
  data_to_dump$END <- as.numeric(data_to_dump$END)

  data_to_dump <- data_to_dump[!is.na(data_to_dump$START),]

  # sort by CHR START and END
  data_to_dump <- data_to_dump[order(data_to_dump$CHR, data_to_dump$START, data_to_dump$END),]

  if (!plyr::empty(data_to_dump)) {

    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " trying to save: ", fileName)

    # save file bed per sample
    # utils::write.table(data_to_dump, file = gzfile(fileName), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    readr::write_tsv(data_to_dump, fileName, col_names = FALSE, na = "NA", progress = FALSE)

    # save file bed per sample temporary to reuse for aggregated bed file
    # filePath <- paste(fileName,"",".temp")
    # sample_names <- rep(sampleName, dim(data_to_dump)[1])
    # data_to_dump <- data.frame(data_to_dump, sample_names)
    # colnames(data_to_dump) <- multipleFileColNames
    #
    # utils::write.table(data_to_dump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " io_dump_sample_as_bed_file: ", fileName)
}


