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
  mat <- as.matrix(signal_data)

  if (ssEnv$inpute == "median")
  {
    # Vectorized: matrixStats::rowMedians is C-level, ~10x faster than apply
    row_med <- matrixStats::rowMedians(mat, na.rm = TRUE)
    na_idx  <- which(is.na(mat), arr.ind = TRUE)
    mat[na_idx] <- row_med[na_idx[, 1L]]
    signal_data <- as.data.frame(mat)
  }
  else if (ssEnv$inpute == "mean")
  {
    row_mn <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    na_idx <- which(is.na(mat), arr.ind = TRUE)
    mat[na_idx] <- row_mn[na_idx[, 1L]]
    signal_data <- as.data.frame(mat)
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
    imputed_matrix <- KMEANS.KNN::kmeans_knn(mat, centers = centers, k = k)
    signal_data <- as.data.frame(imputed_matrix)
  }
  else if (ssEnv$inpute == "none")
  {
    # Vectorized: rowAnys(is.na) avoids creating full logical matrix copy
    has_na <- matrixStats::rowAnyNAs(mat)
    signal_data <- signal_data[!has_na, ]
  }
  else
  {
    log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") ," Invalid inpute value. Please provide one of the following: median, mean, knn.")
    stop()
  }

  rm(mat)
  gc()

  # drop rows with all NA (vectorized)
  all_na <- matrixStats::rowAlls(is.na(as.matrix(signal_data)))
  signal_data <- signal_data[!all_na, ]
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
