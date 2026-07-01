#' sort the dataframe using CHR and START sorting column first for CHR and after
#' for START
#'
#' @param dataframe dataframe to be sorted
#'
#' @return sorted dataframe

#'
anno_sort_by_chr_and_start <- function(dataframe) {

  if (nrow(dataframe) == 0) {
    return(dataframe)
  }

  # TODO: verify CHR and START are column of the data frame
  dataframe$CHR <- as.factor(dataframe$CHR)
  dataframe$CHR <- droplevels(dataframe$CHR)

  dataframe$CHR <- factor(dataframe$CHR, levels = gtools::mixedsort(attributes(dataframe$CHR)$levels))
  dataframe <- dataframe[with(dataframe, order(CHR, START)), ]
  return(dataframe)
}
