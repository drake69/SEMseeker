# Internal utilities for accessing Illumina methylation annotation packages.
#
# The three IlluminaHumanMethylation*anno packages store their data as S4
# objects of class "IlluminaHumanMethylationAnnotation" (defined in minfi).
# Each object has a @data slot containing named data.frames:
#   - Locations    : chr, pos, strand  (rownames = probe IDs)
#   - Islands.UCSC : Islands_Name, Relation_to_Island
#   - Other        : UCSC_RefGene_Name, UCSC_RefGene_Group, ...
#
# SEMseeker accesses these slots directly — minfi is not a dependency.

# Map from SEMseeker technology key to Bioconductor annotation package name
.ANNO_PKGS <- c(
  K850 = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  K450 = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  K27  = "IlluminaHumanMethylation27kanno.ilmn12.hg19"
)

#' Load an Illumina annotation package S4 object
#'
#' The main data object in IlluminaHumanMethylation*anno packages is stored as
#' lazy data, not exported into the namespace.  \code{get(pkg, envir =
#' asNamespace(pkg))} therefore fails.  The correct approach is to use
#' \code{data()} to load the object into a temporary environment.
#'
#' @param pkg Character scalar: annotation package name.
#' @return The \code{IlluminaMethylationAnnotation} S4 object.
#' @keywords internal
.anno_pkg_load <- function(pkg) {
  env <- new.env(parent = emptyenv())
  utils::data(list = pkg, package = pkg, envir = env)
  get(pkg, envir = env)
}

#' Get probe IDs from an Illumina annotation package
#'
#' Extracts the complete vector of probe identifiers (e.g. \code{cg00000029})
#' from an annotation package.
#'
#' @param pkg Character scalar: annotation package name.
#' @return Character vector of probe IDs (rownames of the Locations table).
#' @keywords internal
.anno_pkg_probe_ids <- function(pkg) {
  obj  <- .anno_pkg_load(pkg)
  locs <- slot(obj, "data")$Locations
  rownames(locs)
}

#' Build a data.frame from an Illumina annotation package
#'
#' Combines the Locations, Islands.UCSC, and Other tables from the annotation
#' package S4 object into a single data.frame with one row per probe.
#' Does not require minfi.
#'
#' @param pkg Character scalar: annotation package name.
#' @return A \code{data.frame} with probe IDs as rownames and columns:
#'   \code{chr}, \code{pos}, \code{strand}, \code{Islands_Name},
#'   \code{Relation_to_Island}, \code{UCSC_RefGene_Name},
#'   \code{UCSC_RefGene_Group}, and further columns from the Other table.
#' @keywords internal
.anno_pkg_to_df <- function(pkg) {

  obj       <- .anno_pkg_load(pkg)
  data_list <- slot(obj, "data")

  # Core genomic positions
  locs <- as.data.frame(data_list$Locations, stringsAsFactors = FALSE)

  # CpG island context
  if ("Islands.UCSC" %in% names(data_list)) {
    islands <- as.data.frame(data_list[["Islands.UCSC"]],
                             stringsAsFactors = FALSE)
  } else {
    islands <- data.frame(
      Islands_Name       = NA_character_,
      Relation_to_Island = NA_character_,
      row.names          = rownames(locs),
      stringsAsFactors   = FALSE
    )
  }

  # Gene body and other annotations
  if ("Other" %in% names(data_list)) {
    other <- as.data.frame(data_list$Other, stringsAsFactors = FALSE)
  } else {
    other <- data.frame(row.names = rownames(locs))
  }

  # Combine — all tables share the same rownames (probe IDs)
  df <- cbind(locs, islands, other)
  df
}

#' Detect Illumina array technology by probe-ID overlap
#'
#' Queries each installed annotation package and counts how many probe IDs
#' from \code{probe_ids} are present in each array's probe list.  Returns the
#' technology key (\code{"K27"}, \code{"K450"}, or \code{"K850"}) with the
#' highest overlap count.
#'
#' Returns \code{""} if none of the annotation packages are installed.
#'
#' @param probe_ids Character vector of probe identifiers from the signal matrix.
#' @return Named integer vector of overlap counts, or \code{""} if no packages
#'   are available.
#' @keywords internal
.detect_tech_from_anno <- function(probe_ids) {

  counts <- integer(0)

  for (tech in names(.ANNO_PKGS)) {
    pkg <- .ANNO_PKGS[[tech]]
    if (requireNamespace(pkg, quietly = TRUE)) {
      n <- sum(probe_ids %in% .anno_pkg_probe_ids(pkg))
      counts[[tech]] <- n
    }
  }

  if (length(counts) == 0 || max(counts) == 0)
    return("")

  names(which.max(counts))
}
