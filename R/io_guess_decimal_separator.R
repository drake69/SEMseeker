io_guess_decimal_separator <- function(file) {
  sample_data <- readLines(file, n = 5)
  comma_count <- sum(grepl(",", sample_data))
  period_count <- sum(grepl("\\.", sample_data))
  if (comma_count > period_count) {
    return(",")
  } else {
    return(".")
  }
}
