
position_pivot_to_probe <- function(signal_data)
{
  ssEnv          <- get_session_info()
  probe_features <- probe_annotation_build(ssEnv$tech)[, c("PROBE", "CHR", "START", "END", ssEnv$tech)]
  # transform signal_data into parquet
  signal_data_temp <- polars::as_polars_df(signal_data)
  signal_data_temp <- signal_data_temp$lazy()
  cols_to_remove <- colnames(probe_features)
  cols_to_remove <- cols_to_remove[!(cols_to_remove %in% c("PROBE"))]
  probe_features <- polars::as_polars_df(probe_features)
  probe_features <- probe_features$lazy()
  signal_data_temp <- probe_features$join(
    signal_data_temp,
    on = c("CHR", "START", "END"),
    how = "inner"
  )
  signal_data_temp <- signal_data_temp$unique()$collect()
  signal_data_temp <- as.data.frame(signal_data_temp)
  # Filter to probes matching the detected technology
  tech_col <- ssEnv$tech
  if (tech_col %in% colnames(signal_data_temp))
    signal_data_temp <- signal_data_temp[
      !is.na(signal_data_temp[[tech_col]]) & signal_data_temp[[tech_col]], ]
  signal_data_temp <- signal_data_temp[, !colnames(signal_data_temp) %in% cols_to_remove]
  # move PROBE as first column
  signal_data_temp <- signal_data_temp[,c("PROBE",colnames(signal_data_temp)[!(colnames(signal_data_temp)=="PROBE")])]
  signal_data <- as.data.frame(signal_data_temp)
  rownames(signal_data) <- signal_data$PROBE
  # remove PROBE column
  signal_data$PROBE <- NULL
  return(signal_data)
}
