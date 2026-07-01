#' @title create an annotated file for each marker, figure, area and subarea, each file has all the sample_groups used to calculate epimutation
#'
#' @return nothing
#' @importFrom doRNG %dorng%
anno_annotate_position_pivots <- function ()
{
  start_time <- Sys.time()
  ssEnv <- get_session_info()
  # area and subarea are defined using the filename
  localKeys <-ssEnv$keys_areas_subareas_markers_figures

  # remove POSITION area
  localKeys <- localKeys[localKeys$AREA != "POSITION",]
  # localKeys <- localKeys[localKeys$MARKER != "SIGNAL",]

  if (nrow(localKeys) == 0)
    return()

  # Short-circuit: if every dest pivot already exists on disk, there is
  # nothing to annotate. Avoids the unconditional anno_probe_features_get()
  # load of the Illumina manifest (~10-30s) and the spurious
  # "Annotating genomic area" log line in resume scenarios where no
  # actual annotation work is needed.
  all_dest_exist <- all(vapply(seq_len(nrow(localKeys)), function(i) {
    file.exists(pivot_file_name_parquet(
      localKeys[i, "MARKER"], localKeys[i, "FIGURE"],
      localKeys[i, "AREA"],   localKeys[i, "SUBAREA"]))
  }, logical(1)))
  if (all_dest_exist) {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
      " Annotation skipped: all destination pivots already exist.")
    return()
  }

  log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Annotating genomic area.")

  progress_bar <- ""
  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = seq_len(nrow(localKeys)))

  variables_to_export <- c("ssEnv", "dir_check_and_create", "subarea",
    "progress_bar","progression_index", "progression", "progressor_uuid",
    "owner_session_uuid", "trace","anno_probe_features_get", "localKeys",
    "file_path_build","%>%","get_session_info","log_event")

  # doesn't work with parallel, tests throws error
  for(i in seq_len(nrow(localKeys)))
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
      probe_features <- anno_probe_features_get(area_subarea)
      # NB (AI-027/AI-030): anno_create_position_pivots + stream_merge_bed strip the
      # "chr" prefix from CHR for internal consistency, so the source pivot's
      # CHR is "1", "X" … without prefix. anno_probe_features_get() may return CHR
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

        # AI-027: read via unified dispatcher. CASE 1 (cached parquet) is
        # the normal path here; CASE 2 (streaming merge from bed/bedgraph)
        # makes this resilient to a missing materialised pivot when raw
        # per-sample files still exist.
        pivot <- read_pivot(marker, figure, "POSITION", "WHOLE")
        pivot <- pivot$with_columns(
          polars::pl$col("START")$cast(polars::pl$Int32),
          polars::pl$col("END")$cast(polars::pl$Int32),
          polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", "")
        )
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

        # AI-050: Bioconductor anno-packages assign some probes to multiple
        # genes (intergenic overlaps, antisense, etc), producing composite
        # AREA strings like "NUDT6;SPATA5". Treating the composite as a
        # single gene was a regression that (a) caused apply_stat_model to
        # fail parsing (PVALUE=NA), (b) inflated false positives downstream
        # because a single p-value got smeared across N enrichment hits.
        # Fix: split on ";", explode to N rows, strip whitespace — each
        # multi-mapped probe now contributes separately to every gene's
        # burden, and the group_by below produces clean mono-gene rows.
        #
        # AI-061+ (2026-06-09): extended to also handle "," and "/" as
        # multi-gene separators. Bioconductor manifests use ";" by default
        # but external/custom annotations occasionally use commas, and
        # paralog-cluster compact notations like "HBA1/HBA2" (alpha
        # hemoglobin twin loci) and "HLA-A/B/C" (MHC class I, 3 distinct
        # genes on 6p21.33) are biologically multi-gene — splitting them
        # yields one row per HGNC symbol, which is the correct unit for
        # downstream gene-burden aggregation.
        #
        # AI-107 (2026-06-09): "/" needs SMART splitting with prefix
        # recovery so "HLA-A/B/C" becomes ("HLA-A","HLA-B","HLA-C") rather
        # than ("HLA-A","B","C"). A naive replace_all("/", ";") loses the
        # shared prefix. The expansion is done R-side on the DISTINCT
        # slash-bearing names only (typically a few hundred genes per
        # annotation), then joined back into the lazy pivot before the
        # standard ";" split + explode. ","-bearing names keep their
        # plain semicolon normalisation since each comma-token is already
        # a complete HGNC symbol.
        slashed_areas <- as.data.frame(
          pivot$select("AREA")$
            filter(polars::pl$col("AREA")$str$contains("/", literal = TRUE))$
            unique()$
            collect()
        )$AREA
        if (length(slashed_areas) > 0L) {
          expansions <- vapply(
            slashed_areas,
            function(s) paste(.anno_smart_split_area_name(s), collapse = ";"),
            character(1)
          )
          mapping_lf <- polars::as_polars_df(data.frame(
            AREA      = slashed_areas,
            AREA_NEW  = expansions,
            stringsAsFactors = FALSE
          ))$lazy()
          pivot <- pivot$join(
            mapping_lf, on = "AREA", how = "left"
          )$with_columns(
            polars::pl$when(polars::pl$col("AREA_NEW")$is_not_null())$
              then(polars::pl$col("AREA_NEW"))$
              otherwise(polars::pl$col("AREA"))$alias("AREA")
          )$drop("AREA_NEW")
        }
        pivot <- pivot$with_columns(
          polars::pl$col("AREA")$str$
            replace_all(",", ";")$str$
            split(";")
        )$explode("AREA")
        pivot <- pivot$with_columns(
          polars::pl$col("AREA")$str$strip_chars()
        )
        # Drop rows where AREA became empty after split/strip (e.g. trailing
        # semicolons from malformed annotations).
        pivot <- pivot$filter(polars::pl$col("AREA")$str$len_chars() > 0)

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

