
position_pivot_to_probe <- function(signal_data)
{
  ssEnv          <- get_session_info()
  tech_col       <- ssEnv$tech
  probe_features <- probe_annotation_build(ssEnv$tech)[, c("PROBE", "CHR", "START", "END", tech_col)]
  probe_features_lf <- polars::as_polars_df(probe_features)$lazy()

  # Accept polars LazyFrame, polars DataFrame, or R data.frame/matrix.
  # Staying lazy here avoids materialising the full SIGNAL pivot (e.g. 367k
  # probes × 4k samples ≈ 12 GB) into R memory before the join.
  signal_lazy <-
    if (inherits(signal_data, "polars_lazy_frame"))      signal_data
    else if (inherits(signal_data, "polars_data_frame")) signal_data$lazy()
    else                                                  polars::as_polars_df(signal_data)$lazy()

  # Single lazy chain: join on coords, dedupe, filter probes matching the
  # detected technology, drop annotation columns. ONE collect() at the end.
  result <- probe_features_lf$
    join(signal_lazy, on = c("CHR", "START", "END"), how = "inner")$
    unique()$
    filter(polars::pl$col(tech_col)$is_not_null() & polars::pl$col(tech_col))$
    drop(c("CHR", "START", "END", tech_col))$
    collect()

  signal_data <- as.data.frame(result)
  rownames(signal_data) <- signal_data$PROBE
  signal_data$PROBE <- NULL
  signal_data
}
