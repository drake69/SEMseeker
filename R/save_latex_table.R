save_latex_table <- function(data, dest_file_name_csv, caption)
{
  # remove underscore from everywhere
  colnames(data) <- gsub("_", " ", colnames(data))
  rownames(data) <- gsub("_", " ", rownames(data))
  # remove underscore from everywhere
  data[] <- lapply(data, function(x) gsub("_", " ", x))

  # Convert the dataframe to an xtable object
  table <- xtable::xtable(data)

  # Custom print function to add bold for the first row and add borders
  print.xtable.custom <- function(x, ...) {
    cat("\\begin{table}[H]\n\\centering\n\\scriptsize\n")
    xtable::print.xtable(x,
      include.rownames = FALSE,
      hline.after = c(-1, 0, 1:nrow(x)), # Add horizontal lines after each row
      sanitize.colnames.function = function(y) paste0('\\textbf{', y, '}'),
      booktabs = FALSE,
      floating = FALSE)
    cat(paste0("\\caption{", caption,"}\n"), sep="")
    # create random string
    random_string <- paste(sample(LETTERS, 10, replace = TRUE), collapse = "")
    cat(paste0("\\label{tab:my-table-",random_string,"}\n"), sep="")
    cat("\\end{table}\n")
  }
  # Redirect output to a file
  dest_file_name_tex <- gsub(".csv$", ".tex", dest_file_name_csv)
  sink(dest_file_name_tex)
  print.xtable.custom(table)
  sink() # Res
}
