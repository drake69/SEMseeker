#' Retrieve probe feature annotations for a given genomic area
#'
#' Returns a \code{data.frame} of CpG probe coordinates and feature annotations
#' for the requested area/subarea combination.  Annotations are built from
#' Bioconductor array annotation packages (see \code{\link{probe_annotation_build}})
#' and cached in the session environment, so the package is loaded only once
#' per session.
#'
#' Probes on sex chromosomes are removed when \code{ssEnv$sex_chromosome_remove}
#' is \code{TRUE}.
#'
#' @param area_subarea Character scalar: area and subarea joined by an
#'   underscore (e.g. \code{"GENE_BODY"}, \code{"CHR_WHOLE"},
#'   \code{"ISLAND_N_SHORE"}).  If no underscore is present, \code{"_WHOLE"}
#'   is appended automatically.
#'
#' @return A \code{data.frame} with columns \code{PROBE}, \code{CHR},
#'   \code{START}, \code{END}, and the requested feature column, filtered to
#'   probes matching the current array technology stored in the session
#'   environment.
#'
probe_features_get <- function(area_subarea) {

  ssEnv <- get_session_info()

  # Detect technology if not yet defined
  if (is.null(ssEnv$tech) || ssEnv$tech == "") {
    pivot_file_name <- pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")
    signal_data_pl  <- polars::pl$read_parquet(pivot_file_name)
    signal_data_r   <- as.data.frame(signal_data_pl)
    if ("AREA" %in% colnames(signal_data_r)) {
      rownames(signal_data_r) <- signal_data_r$AREA
      signal_data_r$AREA <- NULL
    }
    ssEnv <- get_meth_tech(signal_data_r)
    log_event("WARNING: probe_features_get() called before technology was defined.")
    if (is.null(ssEnv$tech) || ssEnv$tech == "") {
      log_event("ERROR: could not determine array technology.")
      stop("Could not determine array technology.")
    }
  }

  # Build probe annotation from the Bioconductor array annotation package.
  pkg <- .ANNO_PKGS[[ssEnv$tech]]
  if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) {
    stop("Annotation package '", pkg, "' is not installed. ",
         "Install it with: BiocManager::install('", pkg, "')")
  }
  probe_features <- probe_annotation_build(ssEnv$tech)

  if (!grepl("_", area_subarea))
    area_subarea <- paste0(area_subarea, "_WHOLE")

  # Keep only probes matching the current technology
  probe_features <- probe_features[
    !is.na(probe_features[[ssEnv$tech]]) & probe_features[[ssEnv$tech]], ]
  probe_features$END <- probe_features$START

  if ((grepl("CHR", area_subarea) || grepl("PROBE", area_subarea)) &&
      !grepl("CHR_CYTOBAND", area_subarea)) {
    probe_features <- dplyr::distinct(
      probe_features[, c(ssEnv$tech, "PROBE", "CHR", "START", "END")])
  } else {
    cols_needed <- c(ssEnv$tech, "PROBE", "CHR", "START", "END", area_subarea)
    cols_needed <- intersect(cols_needed, colnames(probe_features))
    probe_features <- dplyr::distinct(probe_features[, cols_needed, drop = FALSE])
  }

  # Add convenience whole-chromosome / probe columns used downstream
  if (grepl("CHR", area_subarea) && !grepl("CHR_CYTOBAND", area_subarea))
    probe_features$CHR_WHOLE <- paste0("chr", probe_features$CHR)

  if (grepl("PROBE", area_subarea))
    probe_features$PROBE_WHOLE <- probe_features$PROBE

  # Drop the technology flag column — not needed downstream
  probe_features <- probe_features[
    , -which(colnames(probe_features) %in% ssEnv$tech), drop = FALSE]

  # Remove sex chromosomes if requested
  if (isTRUE(ssEnv$sex_chromosome_remove))
    probe_features <- probe_features[
      !(probe_features$CHR %in% c("X", "Y")), ]

  return(probe_features)
}
