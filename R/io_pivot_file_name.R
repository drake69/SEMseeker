#' Compose the on-disk name of a pivot (internal)
#'
#' AI-248. **The** single compositor of pivot file names — the file-name twin of
#' [io_feature_colname()], which composes column names. Everything that reads or
#' writes a pivot goes through here, so a pivot's identity is its name.
#'
#' \preformatted{
#' <MARKER>_<FIGURE>_<AREA>_<SUBAREA>[_<AGGREGATION>]_<build>.parquet
#'
#' MUTATIONS_HYPER_CHR_CYTOBAND_SUM_hg19.parquet
#' SIGNAL_BETA_CHR_CYTOBAND_MEAN_hg19.parquet
#' SIGNAL_MVALUE_CHR_CYTOBAND_MEAN_hg19.parquet
#' }
#'
#' The aggregation belongs in the name because before AI-248 the name did **not**
#' say which operator had produced the file: an existing pivot was reused on
#' trust. With one operator per marker that trust was justified; with several it
#' would be the worst kind of bug — silent, because the file is there and the run
#' does not complain while a mean is consumed as if it were a sum. With the
#' aggregation in the key, an existing `MEAN` pivot is a `MEAN` and need not be
#' recomputed, and a missing `MEDIAN` is visible.
#'
#' The same reasoning applies to the `SIGNAL` figure carrying the scale: today a
#' beta run and an M-value run write the very same file name and overwrite each
#' other in the same folder.
#'
#' @param marker,figure,area,subarea the four keys of the pivot.
#' @param aggregation the AGGREGATION axis; `NULL` omits it, which is how the
#'   POSITION-level pivots (one row per position, nothing aggregated yet) are
#'   named.
#' @param add_gz passed through for the CSV flavour.
#' @keywords internal
#' @noRd
io_pivot_file_name <- function(marker, figure, area, subarea, add_gz = TRUE,
                               aggregation = NULL)
{
  .io_pivot_path(marker, figure, area, subarea, aggregation, ".csv", add_gz)
}

#' @rdname io_pivot_file_name
#' @keywords internal
#' @noRd
io_pivot_file_name_parquet <- function(marker, figure, area, subarea,
                                       aggregation = NULL)
{
  .io_pivot_path(marker, figure, area, subarea, aggregation, ".parquet", FALSE)
}

#' @keywords internal
#' @noRd
.io_pivot_path <- function(marker, figure, area, subarea, aggregation,
                           extension, add_gz)
{
  ssEnv <- core_get_session_info()
  reportFolder <- io_dir_check_and_create(ssEnv$result_folderData, "Pivots")
  pivot_subfolder <- io_dir_check_and_create(reportFolder, marker)
  # C-06: append genome_build suffix (e.g. "_hg19") for belt-and-suspenders provenance
  genome_suffix <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    ssEnv$genome_build else "hg19"

  parts <- c(marker, figure, area, subarea)
  if (!is.null(aggregation) && nzchar(as.character(aggregation)))
    parts <- c(parts, as.character(aggregation))
  parts <- c(parts, genome_suffix)

  pivot_base <- paste(parts, collapse = "_")
  io_file_path_build(baseFolder = pivot_subfolder, detailsFilename = pivot_base,
                     extension = extension, add_gz = add_gz)
}
