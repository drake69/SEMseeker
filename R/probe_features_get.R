#' Retrieve probe feature annotations for a given genomic area
#'
#' Returns a \code{data.frame} of CpG probe coordinates and feature annotations
#' for the requested area/subarea combination.
#'
#' For \strong{Illumina} data (K850/K450/K27), annotations are built from
#' Bioconductor array annotation packages (see \code{\link{probe_annotation_build}})
#' and cached in the session environment.
#'
#' For \strong{WGBS} and \strong{LONGREAD} data, coordinates are read directly
#' from the saved POSITION pivot parquet (no Bioconductor annotation needed).
#' Only coordinate-based areas (\code{CHR_WHOLE}, \code{PROBE_WHOLE}) are
#' supported in this release; gene-body and island areas require
#' \code{area_granges_build()} (backlog C-04).
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
#'   \code{START}, \code{END}, and the requested feature column.
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

  if (!grepl("_", area_subarea))
    area_subarea <- paste0(area_subarea, "_WHOLE")

  # -----------------------------------------------------------------------
  # WGBS / LONGREAD path — read positions from the saved POSITION pivot.
  # Only coordinate-based areas are supported; gene/island areas need C-04.
  # -----------------------------------------------------------------------
  if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {

    .coord_areas <- c("CHR", "CHR_WHOLE", "PROBE", "PROBE_WHOLE",
                      "POSITION", "POSITION_WHOLE")
    if (!any(sapply(.coord_areas, function(a) grepl(a, area_subarea,
                                                     fixed = TRUE)))) {
      stop(
        "Area '", area_subarea, "' is not yet supported for ", ssEnv$tech, " data.\n",
        "Gene-body, CpG-island, and other semantic areas require area_granges_build() ",
        "(SEMseeker backlog C-04, planned for a future release).\n",
        "Currently supported areas: CHR_WHOLE, PROBE_WHOLE."
      )
    }

    pf_path <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION", "WHOLE")
    if (!file.exists(pf_path))
      stop("POSITION pivot not found. Ensure signal_save() has completed before ",
           "calling probe_features_get() for area '", area_subarea, "'.")

    pos      <- as.data.frame(polars::pl$read_parquet(pf_path))[, c("CHR","START","END")]
    probe_features <- data.frame(
      PROBE = paste0(pos$CHR, "_", pos$START),
      CHR   = pos$CHR,
      START = pos$START,
      END   = pos$END,
      stringsAsFactors = FALSE
    )

    if (grepl("CHR", area_subarea) && !grepl("CHR_CYTOBAND", area_subarea))
      probe_features$CHR_WHOLE <- paste0("chr", probe_features$CHR)

    if (grepl("PROBE", area_subarea))
      probe_features$PROBE_WHOLE <- probe_features$PROBE

    if (isTRUE(ssEnv$sex_chromosome_remove))
      probe_features <- probe_features[
        !(probe_features$CHR %in% c("X", "Y")), ]

    return(probe_features)
  }

  # -----------------------------------------------------------------------
  # Illumina path — Bioconductor annotation packages
  # -----------------------------------------------------------------------
  pkg <- .ANNO_PKGS[[ssEnv$tech]]
  if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) {
    stop("Annotation package '", pkg, "' is not installed. ",
         "Install it with: BiocManager::install('", pkg, "')")
  }
  probe_features <- probe_annotation_build(ssEnv$tech)

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
