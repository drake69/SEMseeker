#' given data and colnames dump as bed file
#'
#' @param data_to_dump data frame to dump into bed file with CHR, START, END
#' @param fileName name of the file to save data in
#'
#' @return nothing
#'

dump_sample_as_bed_file <- function(data_to_dump, fileName) {

  ssEnv <- get_session_info()
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " dump_sample_as_bed_file ssEnv:", length(ssEnv))
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " dump_sample_as_bed_file:", ssEnv$result_folderData)

  if (!plyr::empty(data_to_dump)) {
    data_to_dump[, "CHR"] <- normalize_chr(data_to_dump[, "CHR"], "output")
  }

  # bed coordinate must start from zero!
  data_to_dump$START <- as.numeric(data_to_dump$START)
  data_to_dump$END <- as.numeric(data_to_dump$END)

  data_to_dump <- data_to_dump[!is.na(data_to_dump$START),]

  # sort by CHR START and END
  data_to_dump <- data_to_dump[order(data_to_dump$CHR, data_to_dump$START, data_to_dump$END),]

  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " trying to save: ", fileName)

  # Always write the BED file, even when data_to_dump has 0 rows. Downstream
  # invariants (e.g. test-6-semseeker.R:83) and consumers assume that the
  # per-sample BED for every (marker, figure) is present on disk whenever
  # the MUTATIONS BED for the same (sample, figure) is present, regardless
  # of whether the derived marker happened to have any rows for that sample.
  # readr::write_tsv on an empty data.frame writes a valid empty .gz.
  readr::write_tsv(data_to_dump, fileName, col_names = FALSE, na = "NA", progress = FALSE)
  log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"),  " dump_sample_as_bed_file: ", fileName)
}


