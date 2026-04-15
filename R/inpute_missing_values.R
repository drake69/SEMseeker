inpute_missing_values <- function(signal_data){

  ssEnv <- get_session_info()
  # count number of na
  n_na <- sum(is.na(signal_data))

  if (n_na==0)
    return (signal_data)

  nrow_ex_ante <- nrow(signal_data)

  # count missed per rows
  n_missed_per_row <- rowSums(is.na(signal_data))
  if (any(n_missed_per_row > 0.1 * ncol(signal_data))) {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " There are rows with missing values more than 10%. I will remove them.")
    log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/(nrow(signal_data)*ncol(signal_data)), 2), " % of the data.")
    signal_data <- signal_data[n_missed_per_row  < 0.1 * ncol(signal_data), ]
  }

  n_missed_per_col <- colSums(is.na(signal_data))
  if (any(n_missed_per_col > 0.1 * nrow(signal_data))) {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " There are columns with missing values more than 10%. I will remove them.")
    log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " There are columns with missing values more than 10%. I will remove them.")
    signal_data <- signal_data[, n_missed_per_col  < 0.1 * nrow(signal_data)]
  }

  chunk_size <- 10000  # Define a chunk size
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_along(seq(1, nrow(signal_data), by = chunk_size)))
  else
    progress_bar <- ""

  n_item = nrow(signal_data)*ncol(signal_data)/100
  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/n_item, 2), " % of the data.")
  log_event("JOURNAL: Imputing missing values using ", ssEnv$inpute , " method. Number of missing values: ", n_na, " corresponding to: ", round(n_na/n_item, 2), " % of the data.")
  if (ssEnv$inpute=="median")
  {
    # Assuming signal_data is a matrix or data frame
    for (i in seq(1, nrow(signal_data), by = chunk_size)) {
      # Define the chunk
      chunk_indices <- i:min(i + chunk_size - 1, nrow(signal_data))
      # Process the chunk
      row_medians <- apply(signal_data[chunk_indices, ], 1, stats::median, na.rm = TRUE)
      row_id <- row(signal_data[chunk_indices, ])[is.na(signal_data[chunk_indices, ])]
      signal_data[chunk_indices, ][is.na(signal_data[chunk_indices, ])] <- row_medians[row_id]
      if(ssEnv$showprogress)
        progress_bar(sprintf("Inputed: %s positions", stringr::str_pad(max(chunk_indices), 10, pad = " ")))

    }
  }
  else if (ssEnv$inpute=="mean")
  {
    for (i in seq(1, nrow(signal_data), by = chunk_size)) {
      # Define the chunk
      chunk_indices <- i:min(i + chunk_size - 1, nrow(signal_data))
      # Process the chunk
      row_medians <- apply(signal_data[chunk_indices, ], 1, mean, na.rm = TRUE)
      signal_data[chunk_indices, ][is.na(signal_data[chunk_indices, ])] <-
        row_medians[row(signal_data[chunk_indices, ])[is.na(signal_data[chunk_indices, ])]]
    }
  }
  else if (grepl("knn", ssEnv$inpute))
  {
    if (length(strsplit(ssEnv$inpute, ";")[[1]]) != 3)
    {
      log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide the number of centers and the number of clusters separated by a semicolon.")
      stop()
    }
    centers <- strsplit(ssEnv$inpute, ";")[[1]][2]
    k <- strsplit(ssEnv$inpute, ";")[[1]][3]
    signal_matrix <- as.matrix(signal_data)
    imputed_matrix <- KMEANS.KNN::kmeans_knn(signal_matrix, centers = centers, k = k)
  }
  else if (ssEnv$inpute=="none")
  {
    signal_data <- signal_data[complete.cases(signal_data),]
  }
  else
  {
    log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide one of the following: median, mean, knn.")
    stop()
  }

  gc()

  # drop rows with all NA
  signal_data <- signal_data[!apply(signal_data, 1, function(x) all(is.na(x))), ]
  nrows_ex_post <- nrow(signal_data)
  if (nrows_ex_post < nrow_ex_ante)
  {
    log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Dropping rows with all missing values. Number of rows dropped: ", nrow_ex_ante - nrows_ex_post)
    log_event("JOURNAL: Dropped rows with all missing values. Number of rows dropped: ", nrow_ex_ante - nrows_ex_post)
  }
  n_na <- sum(is.na(signal_data))
  log_event("INFO:", format(Sys.time(), "%a %b %d %X %Y") ," Imputation completed. Number of missing values: ", n_na)
  return(signal_data)
}
