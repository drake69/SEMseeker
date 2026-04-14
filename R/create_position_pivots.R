#' create pivot from samples' single bed for the PROBE (WHOLE)
#'
#' @param population population data frame to use for creating pivot
#' @param keys markers and figure to create pivot of
#'
#' @return A data.frame of genomic position pivots with chromosomal coordinates
#'   and per-sample signal values; results are also written to parquet files on
#'   disk. Returns invisibly.
#'
#' @importFrom doRNG %dorng%
create_position_pivots <- function(population, keys)
{
  ssEnv <- get_session_info()

  # create pivot for basic markers and figure
  ssEnv <- get_session_info()
  selection <-c("MUTATIONS","DELTAR","DELTAS")
  # sort keys by MARKER
  keys <- keys[order(keys$MARKER),]
  keys <- subset(keys, MARKER %in% selection)
  population_result <- data.frame()
  if(nrow(keys)==0)
    return()

  variables_to_export <- c("pivot", "pivot_filename", "temp_pop", "samples", "progress_bar_pop", "progress_bar_key", "progress_bar_pop",
    "sample", "sample_id", "sample_group", "pivot_filename", "marker", "figure", "area", "subarea", "s", "progress_bar_pop")

  # foreach::foreach(k=1:nrow(keys), .export = variables_to_export) %dorng%
  for (k in 1:nrow(keys))
  {
    key <- keys[k,]
    marker <- as.character(key$MARKER)
    figure <- as.character(key$FIGURE)
    area <- as.character("POSITION")
    subarea <- as.character("WHOLE")
    if(is.na(marker) || is.na(figure) || is.na(area) || is.na(subarea))
      next
    pivot_filename <- pivot_file_name_parquet(marker, figure, area, subarea)
    # look if pivot exists
    # if(ssEnv$showprogress)
    #   progress_bar_key(sprintf("Creating pivot of %s %s %s %s", marker, figure, area, subarea))
    # read first row of parquet file
    temp_pop <- population
    if (file.exists(pivot_filename))
    {
      temp <- as.data.frame(polars::pl$scan_parquet(pivot_filename, n_rows=1)$collect())
      samples <- colnames(temp)
      temp_pop <- temp_pop[!(temp_pop$Sample_ID %in% samples),]
    }
    if(ssEnv$showprogress)
      progress_bar_pop <- progressr::progressor(along = 1:(nrow(temp_pop)))
    else
      progress_bar_pop <- ""

    # remove rows with empty Sample_Group
    temp_pop <- temp_pop[!is.na(temp_pop$Sample_Group),]
    if(nrow(temp_pop) == 0)
      next
    for ( s in 1:nrow(temp_pop))
    {

      gc()
      if(ssEnv$showprogress)
        progress_bar_pop(sprintf("Creating pivot of %s %s %s %s", marker, figure, area, subarea))

      sample <- temp_pop[s,]
      sample_id <- as.character(sample$Sample_ID)
      sample_group <- as.character(sample$Sample_Group)
      if(!exists("pivot"))
        if (file.exists(pivot_filename))
          pivot <- polars::pl$scan_parquet(pivot_filename)

      # pivot_colnames <- colnames(pivot)
      # if(sample_id %in% pivot_colnames)
      #   next

      # read bed file
      bed_file <- bed_file_name(sample_id,sample_group,marker,figure)
      if (!file.exists(bed_file))
        next
      else
      {
        # Read the file in lazy mode
        # data <- polars::pl$scan_csv(bed_file, has_header = FALSE, separator = "\t", skip_rows = 0)
        data <- readr::read_tsv(bed_file,
          col_types = readr::cols(
            .default = readr::col_double(),
            X1 = readr::col_character(),
            X2 = readr::col_integer(),
            X3 = readr::col_integer()
          ),
          show_col_types=FALSE, progress=FALSE, skip=0, col_names=FALSE)

        data <- polars::as_polars_df(data)
        data <- data$lazy()

        # 4 columns → Rename normally (use select+alias for cross-version polars compatibility)
        data <- data$select(
          polars::pl$col("X1")$alias("CHR"),
          polars::pl$col("X2")$alias("START"),
          polars::pl$col("X3")$alias("END"),
          polars::pl$col("X4")$alias(sample_id)
        )

        # Apply filters
        data <- data$filter(
          polars::pl$col("CHR")$is_not_null() &
            polars::pl$col("START")$is_not_null() &
            polars::pl$col("END")$is_not_null() &
            polars::pl$col("START") > 0 &
            polars::pl$col("END") > 0 &
            polars::pl$col("START") <= polars::pl$col("END")
        )
      }

      if(!exists("pivot"))
      {
        pivot <- data
        rm(data)
        next
      }
      pivot <- pivot$join(
        data,
        on = c("CHR", "START", "END"),
        how = "full"
      )
      # Example: Assuming pivot is already the result of your previous join
      pivot <- pivot$with_columns(
        polars::pl$when(polars::pl$col("CHR")$is_not_null())$then(polars::pl$col("CHR"))$otherwise(polars::pl$col("CHR_right"))$alias("CHR"),
        polars::pl$when(polars::pl$col("START")$is_not_null())$then(polars::pl$col("START"))$otherwise(polars::pl$col("START_right"))$alias("START"),
        polars::pl$when(polars::pl$col("END")$is_not_null())$then(polars::pl$col("END"))$otherwise(polars::pl$col("END_right"))$alias("END")
      )
      pivot <- pivot$drop(c("CHR_right", "START_right","END_right"))
      # message(pivot$columns)
      if ( s %% 100 == 0 )
      {
        pivot <- pivot$collect()
        pivot <- pivot$sort(c("CHR", "START"), descending = FALSE)
        pivot$write_parquet(pivot_filename)
        pivot <- pivot$lazy()
      }
    }

    if (!exists("pivot"))
      next

    pivot <- pivot$sort(c("CHR", "START"), descending = FALSE)
    pivot <- pivot$collect()
    pivot$write_parquet(pivot_filename)
    pivot_sidecar_write(pivot_filename)   # C-06: write genome_build + tech sidecar
    rm(pivot)
  }
}


