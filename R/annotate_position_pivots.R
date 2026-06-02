#' @title create an annotated file for each marker, figure, area and subarea, each file has all the sample_groups used to calculate epimutation
#'
#' @return nothing
#' @importFrom doRNG %dorng%
annotate_position_pivots <- function ()
{
  start_time <- Sys.time()
  ssEnv <- get_session_info()
  # area and subarea are defined using the filename
  localKeys <-ssEnv$keys_areas_subareas_markers_figures

  # remove POSITION area
  localKeys <- localKeys[localKeys$AREA != "POSITION",]
  # localKeys <- localKeys[localKeys$MARKER != "SIGNAL",]

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating genomic area.")

  progress_bar <- ""
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:(nrow(localKeys)))

  variables_to_export <- c("ssEnv", "dir_check_and_create", "subarea",
    "progress_bar","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","probe_features_get", "localKeys",
    "file_path_build","%>%","get_session_info","log_event")

  # check probe features are avaialable
  probe_features <- probe_features_get("PROBE_WHOLE")

  if(nrow(localKeys)==0)
    return()

  # doesn't work with parallel, tests throws error
  for(i in 1:nrow(localKeys))
    # foreach::foreach(i=1:nrow(localKeys), .export = variables_to_export) %dorng%
  {
    marker <- as.character(localKeys[i,"MARKER"])
    figure <- as.character(localKeys[i,"FIGURE"])
    subarea <- as.character(localKeys[i,"SUBAREA"])
    area <- as.character(localKeys[i,"AREA"])
    # TO DO: remove subarea whole from probe features
    area_subarea <- paste(area,"_", ifelse (subarea=="","WHOLE",subarea) , sep="")
    source_pivot_filename <- pivot_file_name_parquet(marker, figure, "POSITION","WHOLE")
    dest_pivot_filename <- pivot_file_name_parquet(marker, figure, area,subarea)
    # i <- 1
    if (!file.exists(dest_pivot_filename))
    {
      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " File does not exists: ", dest_pivot_filename)
      probe_features <- probe_features_get(area_subarea)
      # NB (AI-027/AI-030): create_position_pivots + stream_merge_bed strip the
      # "chr" prefix from CHR for internal consistency, so the source pivot's
      # CHR is "1", "X" … without prefix. probe_features_get() may return CHR
      # either with or without prefix depending on the manifest table; we
      # normalise both sides to no-prefix before the inner join to avoid
      # 0-row joins.

      # annotate file
      if(file.exists(source_pivot_filename))
      {
        probe_features$CHR <- as.character(probe_features$CHR)
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating, reading pivot.")
        probe_features <- polars::as_polars_df(probe_features)$lazy()
        probe_features <- probe_features$with_columns(polars::pl$col(area_subarea)$alias("AREA"))$drop(area_subarea)
        # colnames(probe_features)
        probe_features <- probe_features$with_columns(
          polars::pl$col("START")$cast(polars::pl$Int32),
          polars::pl$col("END")$cast(polars::pl$Int32),
          polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", "")
        )

        pivot <- polars::pl$scan_parquet(source_pivot_filename)
        pivot <- pivot$with_columns(
          polars::pl$col("START")$cast(polars::pl$Int32),
          polars::pl$col("END")$cast(polars::pl$Int32),
          polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", "")
        )
        # pivot <- polars::pl$read_parquet(source_pivot_filename)
        pivot <- probe_features$join(
          pivot,
          on = c("CHR", "START", "END"),
          how = "inner"
        )

        existing_cols <- names(pivot)
        cols_to_remove <- c("PROBE","CHR","START","END","K27","K450","K850")
        cols_to_remove <- cols_to_remove[cols_to_remove %in% existing_cols]
        pivot <- pivot$drop(cols_to_remove)
        # drop row where AREA is NA
        pivot <- pivot$drop_nulls("AREA")
        pivot <- pivot$sort(c("AREA"), descending = FALSE)$collect()

        if (localKeys[i, "DISCRETE"]) {
          pivot <- pivot$group_by("AREA", .maintain_order=FALSE)$sum()
        } else {
          pivot <- pivot$group_by("AREA", .maintain_order=FALSE)$mean()
        }

        pivot$write_parquet(dest_pivot_filename)

        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating, annotaion executed.")

        if(nrow(pivot)==0)
          ssEnv$key_missed_areas_subareas <- unique(rbind(ssEnv$key_missed_areas_subareas, localKeys[i,c("AREA","SUBAREA")]))
      }
    }


    if(ssEnv$showprogress)
      progress_bar(sprintf("Annotating position pivots."))
  }

  # remove missed keys
  selector <- !((ssEnv$keys_areas_subareas_markers_figures$AREA %in% ssEnv$key_missed_areas_subareas$AREA) & (ssEnv$keys_areas_subareas_markers_figures$SUBAREA %in% ssEnv$key_missed_areas_subareas$SUBAREA))
  ssEnv$keys_areas_subareas_markers_figures  <- ssEnv$keys_areas_subareas_markers_figures[selector,]

  update_session_info(ssEnv)
  end_time <- Sys.time()
  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotation genomic areas file finished in ", difftime(end_time,start_time,units = "mins")," minutes.")
  gc()
  #
}

