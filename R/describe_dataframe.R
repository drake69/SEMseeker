#' Describe a data frame with summary statistics per column
#'
#' Pure helper: no I/O, no side effects. Returns a data frame with one row per
#' column of \code{df} and the following fields: Variable, Class, Missing_Values,
#' Missing_Values_Percent, Unique_Values, Mean, Median, Min, Max.
#'
#' @param df A data frame to describe.
#' @return A data frame with summary statistics.
#'
describe_dataframe <- function(df) {
  data.frame(
    Variable              = names(df),
    Class                 = sapply(df, class),
    Missing_Values        = sapply(df, function(x) sum(is.na(x))),
    Missing_Values_Percent = round(sapply(df, function(x) sum(is.na(x)) / length(x) * 100), 2),
    Unique_Values         = sapply(df, function(x) length(unique(x))),
    Mean   = round(sapply(df, function(x) if (is.numeric(x)) mean(x,   na.rm = TRUE) else NA), 2),
    Median = round(sapply(df, function(x) if (is.numeric(x)) median(x, na.rm = TRUE) else NA), 2),
    Min    = round(sapply(df, function(x) if (is.numeric(x)) min(x,    na.rm = TRUE) else NA), 2),
    Max    = round(sapply(df, function(x) if (is.numeric(x)) max(x,    na.rm = TRUE) else NA), 2)
  )
}
