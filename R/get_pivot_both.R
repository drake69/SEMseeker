get_pivot_both <- function(marker)
{
  ssEnv <- get_session_info()
  is_discrete <- unique(ssEnv$keys_markers_figures_default[ssEnv$keys_markers_figures_default$MARKER==marker,"DISCRETE"])

  # Migrated 2026-06-01 to use the new read_pivot() dispatcher (AI-027).
  # Branch (a): BOTH already materialised (cache).
  pivot_file_name_both <- pivot_file_name_parquet(marker,"BOTH","PROBE","WHOLE")
  if(file.exists(pivot_file_name_both))
    return(as.data.frame(polars::pl$read_parquet(pivot_file_name_both)))

  # Branch (b): build BOTH from HYPER+HYPO. read_pivot() transparently handles
  # both the cached-parquet and the bed-streaming-merge backends; if neither
  # the parquet nor the per-sample bed files exist for a figure, it returns
  # NULL and we treat that figure as empty.
  pivot_hyper <- read_pivot(marker, "HYPER", area = "PROBE", subarea = "WHOLE")
  pivot_hypo  <- read_pivot(marker, "HYPO",  area = "PROBE", subarea = "WHOLE")

  if (is.null(pivot_hyper) && is.null(pivot_hypo))
    return(data.frame())
  # Drop the NULL side(s); polars concat requires non-NULL operands
  parts <- Filter(Negate(is.null), list(pivot_hyper, pivot_hypo))

  # Union (concatenate) the available figures.
  # NOTA polars 1.11 (e già da 1.x): pl$concat() vuole i frame come varargs (...),
  # non come list. Usiamo do.call per espandere la lista in posizionali.
  pivot_both <- do.call(polars::pl$concat, parts)$collect()
  pivot_both <- pivot_both$group_by("AREA", .maintain_order=FALSE)$sum()

  pivot_both$write_parquet(pivot_file_name_both)
  # Sidecar JSON is materialised by ensure_sidecars() at pipeline end.

  # pivot_both <- pivot_both$group_by("AREA", maintain_order=FALSE)$agg(
  #   pl$all()$exclude("AREA")$sum()
  # )
  #
  # pivot_both$sink_parquet(pivot_file_name_both)


  # pivot_hyper <- data.frame()
  # if(file.exists(pivot_file_name_hyper))
  #   pivot_hyper <- polars::pl$read_parquet(pivot_file_name_hyper)$to_data_frame()
  # pivot_hypo <- data.frame()
  # if(file.exists(pivot_file_name_hypo))
  #   pivot_hypo <- polars::pl$read_parquet(pivot_file_name_hypo)$to_data_frame()
  #
  # # get the row with the max count of values
  # count_m <- plyr::rbind.fill(as.data.frame(pivot_hyper), as.data.frame(pivot_hypo))
  # rm(pivot_hypo)
  # rm(pivot_hyper)
  #
  #
  # # sort count_m by AREA
  # count_m <- count_m[order(count_m$AREA),]
  #
  # if (nrow(count_m)!=0)
  #   # Group by AREA and sum all other columns
  #   df_grouped <- df$group_by("AREA")$agg(
  #     pl$all()$exclude("AREA")$sum()
  #   )
  #
  # count_m <- as.data.frame(count_m)
  # arrow::write_parquet(count_m, pivot_file_name_both)
  return(as.data.frame(pivot_both))
}
