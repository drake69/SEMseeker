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
  # WGBS / LONGREAD path
  # Coordinate-only areas (CHR_WHOLE, PROBE_WHOLE) are handled inline.
  # All semantic areas (GENE_*, ISLAND_*, CHR_CYTOBAND, DMR_*) are built
  # via area_granges_build() and CpG assignments resolved with findOverlaps().
  # -----------------------------------------------------------------------
  if (ssEnv$tech %in% c("WGBS", "LONGREAD")) {

    pf_path <- pivot_file_name_parquet("SIGNAL", "MEAN", "POSITION", "WHOLE")
    if (!file.exists(pf_path))
      stop("POSITION pivot not found. Ensure signal_save() completed before ",
           "calling probe_features_get() for area '", area_subarea, "'.")

    pos <- as.data.frame(
      polars::pl$read_parquet(pf_path))[, c("CHR", "START", "END")]
    probe_ids <- paste0(pos$CHR, "_", pos$START)

    # --- Coordinate-only areas (no annotation needed) ---
    if (area_subarea %in% c("CHR_WHOLE", "PROBE_WHOLE",
                             "POSITION_WHOLE", "PROBE")) {
      probe_features <- data.frame(
        PROBE = probe_ids, CHR = pos$CHR, START = pos$START, END = pos$END,
        stringsAsFactors = FALSE)
      if (grepl("CHR", area_subarea))
        probe_features$CHR_WHOLE <- paste0("chr", probe_features$CHR)
      if (grepl("PROBE", area_subarea))
        probe_features$PROBE_WHOLE <- probe_features$PROBE
      if (isTRUE(ssEnv$sex_chromosome_remove))
        probe_features <- probe_features[
          !(probe_features$CHR %in% c("X", "Y")), ]
      return(probe_features)
    }

    # --- Semantic areas via area_granges_build() + findOverlaps() ---
    for (pkg in c("GenomicRanges", "IRanges", "S4Vectors")) {
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' is required for WGBS/LONGREAD area analysis.\n",
             "Install: BiocManager::install(c('GenomicRanges','IRanges','S4Vectors'))")
    }

    # Build GRanges for the CpG positions (seqnames need "chr" prefix for TxDb)
    cpg_gr <- GenomicRanges::GRanges(
      seqnames = paste0("chr", pos$CHR),
      ranges   = IRanges::IRanges(pos$START, pos$START)
    )

    # Build area GRanges (cached after first call per session)
    area_gr <- area_granges_build(area_subarea,
                                  genome_build = ssEnv$genome_build)

    # Assign each CpG to its overlapping area (ignore strand for methylation)
    hits <- GenomicRanges::findOverlaps(cpg_gr, area_gr,
                                        ignore.strand = TRUE)

    q_hits <- S4Vectors::queryHits(hits)
    s_hits <- S4Vectors::subjectHits(hits)

    probe_features <- data.frame(
      PROBE = probe_ids[q_hits],
      CHR   = pos$CHR[q_hits],
      START = pos$START[q_hits],
      END   = pos$END[q_hits],
      stringsAsFactors = FALSE
    )
    probe_features[[area_subarea]] <-
      as.character(GenomicRanges::mcols(area_gr)$label[s_hits])

    if (nrow(probe_features) == 0)
      log_event("WARNING: no CpG positions overlap area '", area_subarea,
                "' for genome_build = '", ssEnv$genome_build, "'.")

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
