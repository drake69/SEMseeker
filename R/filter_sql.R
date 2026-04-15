filter_sql <- function(sql_conditions, data_frame)
{

  # check if sql_conditions is am empty string
  if (length(sql_conditions) == 0 )
    return(data_frame)

  for ( s in seq_along(sql_conditions))
  {
    sql_condition <- as.character(sql_conditions[s])
    sql_condition <- trimws(sql_condition)
    if(!is.na(sql_condition))
      if(!is.null(sql_condition))
        if (sql_condition != "")
        {
          # sql_condition <- toupper(sql_condition)
          sql_condition <- gsub("  ", " ", sql_condition)
          if (grepl("FROM", sql_condition))
          {
            if (!grepl("FROM TABLE", sql_condition))
            {
              log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," the only allowed table placeholder with sub statement is TABLE!")
              stop()
            }
            else
              sql_condition <- gsub("TABLE","data_frame", sql_condition)
          }
          columns <- toupper(colnames(data_frame))
          # drop duplicates columns
          data_frame <- data_frame[, !duplicated(columns), drop = FALSE]
          columns <- colnames(data_frame)
          # Use an explicit environment with a simple table name so sqldf
          # can reliably find the data frame regardless of calling context.
          assign("ssdf_tmp", data_frame, envir = .GlobalEnv)
          sql <- paste("select * from ssdf_tmp where ", sql_condition)
          data_frame <- sqldf::sqldf(sql, envir = .GlobalEnv)
          rm("ssdf_tmp", envir = .GlobalEnv)
          log_event("DEBUG:", format(Sys.time(), "%a %b %d %X %Y"), " Executed sql: " , sql_condition)
        }
  }
  return(data_frame)
}
