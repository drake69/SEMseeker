#' Build the probe annotation table from Bioconductor packages
#'
#' Constructs a \code{data.frame} equivalent to the legacy \code{pp_tot} object
#' by reading coordinates and feature annotations directly from the S4 data
#' slots of the installed Illumina annotation package.  The annotation packages
#' are accessed directly via their S4 data slots.
#' The result is cached inside the session environment
#' (\code{ssEnv$probe_annotation}) so the package is parsed only once per
#' R session.
#'
#' Column mapping from Bioconductor to SEMseeker:
#' \tabular{ll}{
#'   \code{PROBE}             \tab rownames of the annotation Locations table  \cr
#'   \code{CHR}               \tab \code{chr} column (strip "chr" prefix)      \cr
#'   \code{START / END}       \tab \code{pos} column (single-base CpG probes)  \cr
#'   \code{GENE_*}            \tab parsed from \code{UCSC_RefGene_Group}       \cr
#'   \code{ISLAND_*}          \tab recoded from \code{Relation_to_Island}      \cr
#'   \code{CHR_CYTOBAND}      \tab \code{UCSC_CpG_Islands_Name} or similar     \cr
#'   \code{DMR_WHOLE/DMR_DMR} \tab from bundled \code{dmr_annotation} table    \cr
#' }
#'
#' @param tech Character scalar: one of \code{"K27"}, \code{"K450"}, or
#'   \code{"K850"}.
#' @param force Logical: if \code{TRUE}, rebuild the cache even when a cached
#'   version already exists.
#'
#' @return A \code{data.frame} with one row per probe and columns
#'   \code{PROBE}, \code{CHR}, \code{START}, \code{END}, the technology flag,
#'   gene body region columns, CpG island context columns, cytoband, and DMR
#'   annotations.
#'
#' Build the eight GENE_* annotation columns from the UCSC RefGene manifest.
#'
#' Pure (no annotation package). \code{GENE_<region>} lists the unique gene
#' symbols whose RefGene group matches that region; \code{GENE_WHOLE} lists all
#' genes overlapping the probe (whole gene), mirroring \code{ISLAND_WHOLE}.
#'
#' @param group_str Character vector: \code{UCSC_RefGene_Group} (";"-joined).
#' @param name_str Character vector: \code{UCSC_RefGene_Name} (";"-joined).
#' @return Named list of GENE_BODY, GENE_TSS200, GENE_TSS1500, GENE_1STEXON,
#'   GENE_5UTR, GENE_3UTR, GENE_EXONBND, GENE_WHOLE.
#' @keywords internal
#' @noRd
.anno_gene_columns <- function(group_str, name_str) {
  gene_region_map <- c(
    GENE_BODY    = "Body",    GENE_TSS200  = "TSS200",  GENE_TSS1500 = "TSS1500",
    GENE_1STEXON = "1stExon", GENE_5UTR    = "5'UTR",   GENE_3UTR    = "3'UTR",
    GENE_EXONBND = "ExonBnd"
  )
  gene_groups <- strsplit(as.character(group_str), ";", fixed = TRUE)
  gene_names  <- strsplit(as.character(name_str),  ";", fixed = TRUE)

  extract <- function(region) {
    vapply(seq_along(gene_groups), function(i) {
      g <- gene_groups[[i]]; n <- gene_names[[i]]
      hits <- unique(n[g == region & n != ""])
      if (length(hits) == 0L) NA_character_ else paste(hits, collapse = ";")
    }, character(1))
  }

  out <- lapply(gene_region_map, extract)  # names() = GENE_BODY, GENE_TSS200, ...
  out$GENE_WHOLE <- vapply(gene_names, function(genes) {
    hits <- unique(genes[genes != "" & !is.na(genes)])
    if (length(hits) == 0L) NA_character_ else paste(hits, collapse = ";")
  }, character(1))
  out
}

#' Assign each probe its cytoband by range overlap against \code{cytoband_hg19}.
#'
#' Pure: one \code{findInterval} per chromosome (O(n log m)), no package access.
#'
#' @param chr Character vector: chromosome WITHOUT the "chr" prefix.
#' @param start Integer vector: CpG position.
#' @param cytoband Cytoband table (CHR/START/END/CYTOBAND); defaults to the
#'   bundled \code{cytoband_hg19} (injectable for testing).
#' @return Named list with a single \code{CHR_CYTOBAND} character vector.
#' @keywords internal
#' @noRd
.anno_chr_columns <- function(chr, start, cytoband = NULL) {
  if (is.null(cytoband)) cytoband <- SEMseeker::cytoband_hg19
  cb <- cytoband[!is.na(cytoband$CHR) & cytoband$CHR != "", , drop = FALSE]
  chr_vec   <- as.character(chr)
  start_vec <- as.integer(start)
  cytoband_vec <- rep(NA_character_, length(chr_vec))

  for (chr_val in unique(chr_vec[!is.na(chr_vec)])) {
    cb_chr <- cb[cb$CHR == chr_val, , drop = FALSE]
    if (nrow(cb_chr) == 0L) next
    cb_chr <- cb_chr[order(cb_chr$START), ]
    idx <- which(chr_vec == chr_val)
    band_idx <- findInterval(start_vec[idx], cb_chr$START)
    valid <- band_idx > 0L & start_vec[idx] <= cb_chr$END[band_idx]
    cytoband_vec[idx[valid]] <- cb_chr$CYTOBAND[band_idx[valid]]
  }
  list(CHR_CYTOBAND = cytoband_vec)
}

#' Build the Illumina probe annotation table
#'
#' Internal helper. Assembles the per-probe annotation (genomic position,
#' cytoband and DMR/area membership) for an Illumina methylation array
#' platform, joining the bundled \code{\link{cytoband_hg19}} and
#' \code{\link{dmr_annotation}} reference data.
#'
#' @param tech Illumina platform identifier (e.g. "EPIC", "450k", "27k").
#' @param force Logical; rebuild even when a cached annotation is available.
#' @return A data frame of per-probe annotation columns.
#' @keywords internal
anno_probe_annotation_build <- function(tech, force = FALSE) {

  ssEnv <- get_session_info()
 
  # Return cached version unless forced
  if (!force &&
      !is.null(ssEnv$probe_annotation) &&
      identical(ssEnv$probe_annotation_tech, tech)) {
    return(ssEnv$probe_annotation)
  }

  pkg <- .ANNO_PKGS[[tech]]
  if (is.null(pkg))
    stop("Unknown technology: '", tech, "'. Must be one of: ",
         paste(names(.ANNO_PKGS), collapse = ", "))

  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Annotation package '", pkg, "' is required but not installed.\n",
         "Install with: BiocManager::install('", pkg, "')")

  log_event("INFO: building probe annotation from '", pkg, "'...")

  # ---- Read raw data from annotation package S4 slots ----
  anno_df <- .anno_pkg_to_df(pkg)

  # ---- PROBE, CHR, START, END ----
  anno_df$PROBE <- rownames(anno_df)
  anno_df$CHR   <- sub("^chr", "", as.character(anno_df$chr))
  anno_df$START <- as.integer(anno_df$pos)
  anno_df$END   <- as.integer(anno_df$pos)

  # ---- Technology flag ----
  anno_df[[tech]] <- TRUE

  # ---- Semantic area columns (one row per probe) ----
  # GENE / ISLAND / CHR are 1:1 mappings: each pure helper returns a NAMED LIST
  # of columns for ALL probes. They are independent column-groups, NOT a
  # mutually-exclusive dispatch — every probe gets its gene context AND its
  # island context AND its cytoband. The helpers are pure (no annotation-package
  # access), so each area's recoding is unit-tested without an Illumina package.
  col_groups <- c(
    .anno_gene_columns(anno_df$UCSC_RefGene_Group, anno_df$UCSC_RefGene_Name),
    .anno_island_columns(anno_df$Relation_to_Island, anno_df$Islands_Name,
                    anno_df$CHR, anno_df$START),
    .anno_chr_columns(anno_df$CHR, anno_df$START)
  )
  for (nm in names(col_groups)) anno_df[[nm]] <- col_groups[[nm]]

  # ---- DMR columns (1:many membership — NOT a 1:1 column) ----
  # A probe can belong to several DMRs, so this is a row-EXPANDING join, not a
  # per-probe column like GENE/ISLAND/CHR. The duplication is intentional and
  # required: anno_probe_features_get() selects [tech, PROBE, CHR, START, END,
  # area_subarea] and dplyr::distinct()s — for DMR_* this preserves every
  # membership, while for the other areas the duplicate rows collapse back.
  dmr <- SEMseeker::dmr_annotation
  anno_df <- merge(anno_df, dmr, by = "PROBE", all.x = TRUE)

  # ---- Select final columns ----
  keep <- c(
    "PROBE", "CHR", "START", "END", tech,
    "GENE_BODY", "GENE_TSS200", "GENE_TSS1500", "GENE_1STEXON",
    "GENE_5UTR", "GENE_3UTR", "GENE_EXONBND", "GENE_WHOLE",
    "ISLAND_WHOLE", "ISLAND_ISLAND",
    "ISLAND_N_SHORE", "ISLAND_S_SHORE",
    "ISLAND_N_SHELF", "ISLAND_S_SHELF", "ISLAND_OPENSEA",
    "CHR_CYTOBAND", "DMR_WHOLE", "DMR_DMR"
  )
  anno_df <- anno_df[, intersect(keep, colnames(anno_df)), drop = FALSE]

  # ---- Cache and return ----
  ssEnv$probe_annotation      <- anno_df
  ssEnv$probe_annotation_tech <- tech
  update_session_info(ssEnv)

  log_event("INFO: probe annotation built — ", nrow(anno_df), " probes, tech = ", tech)
  return(anno_df)
}
