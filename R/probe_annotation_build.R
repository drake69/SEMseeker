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
probe_annotation_build <- function(tech, force = FALSE) {

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

  # ---- GENE columns (vectorized via data.table for speed at 485k+ probes) ----
  group_str <- as.character(anno_df$UCSC_RefGene_Group)
  name_str  <- as.character(anno_df$UCSC_RefGene_Name)

  gene_region_map <- c(
    GENE_BODY    = "Body",
    GENE_TSS200  = "TSS200",
    GENE_TSS1500 = "TSS1500",
    GENE_1STEXON = "1stExon",
    GENE_5UTR    = "5'UTR",
    GENE_3UTR    = "3'UTR",
    GENE_EXONBND = "ExonBnd"
  )

  # Pre-split once (not per-column)
  gene_groups <- strsplit(group_str, ";", fixed = TRUE)
  gene_names  <- strsplit(name_str,  ";", fixed = TRUE)

  # Helper: extract unique gene names matching a region, vectorized via vapply
  .extract_genes_for_region <- function(region, groups_list, names_list) {
    vapply(seq_along(groups_list), function(i) {
      g <- groups_list[[i]]
      n <- names_list[[i]]
      hits <- unique(n[g == region & n != ""])
      if (length(hits) == 0L) NA_character_ else paste(hits, collapse = ";")
    }, character(1))
  }

  for (col in names(gene_region_map)) {
    anno_df[[col]] <- .extract_genes_for_region(
      gene_region_map[[col]], gene_groups, gene_names)
  }

  anno_df$GENE_WHOLE <- vapply(gene_names, function(genes) {
    hits <- unique(genes[genes != "" & !is.na(genes)])
    if (length(hits) == 0L) NA_character_ else paste(hits, collapse = ";")
  }, character(1))

  # ---- ISLAND columns ----
  island_rel <- as.character(anno_df$Relation_to_Island)
  island_name <- as.character(anno_df$Islands_Name)

  anno_df$ISLAND_WHOLE   <- ifelse(island_rel == "Island",  island_name, NA_character_)
  anno_df$ISLAND_N_SHORE <- ifelse(island_rel == "N_Shore", island_name, NA_character_)
  anno_df$ISLAND_S_SHORE <- ifelse(island_rel == "S_Shore", island_name, NA_character_)
  anno_df$ISLAND_N_SHELF <- ifelse(island_rel == "N_Shelf", island_name, NA_character_)
  anno_df$ISLAND_S_SHELF <- ifelse(island_rel == "S_Shelf", island_name, NA_character_)

  # ---- CHR_CYTOBAND ----
  # Assigned by range overlap against the bundled cytoband_hg19 table.
  # Vectorized: one findInterval per chromosome instead of nested loops.
  cb        <- SEMseeker::cytoband_hg19
  cb        <- cb[!is.na(cb$CHR) & cb$CHR != "", ]
  chr_vec   <- anno_df$CHR
  start_vec <- anno_df$START
  cytoband_vec <- rep(NA_character_, nrow(anno_df))

  for (chr_val in unique(chr_vec[!is.na(chr_vec)])) {
    cb_chr <- cb[cb$CHR == chr_val, , drop = FALSE]
    if (nrow(cb_chr) == 0L) next
    cb_chr <- cb_chr[order(cb_chr$START), ]
    idx_anno <- which(chr_vec == chr_val)
    # findInterval: O(n log m) instead of O(n × m)
    band_idx <- findInterval(start_vec[idx_anno], cb_chr$START)
    valid <- band_idx > 0L & start_vec[idx_anno] <= cb_chr$END[band_idx]
    cytoband_vec[idx_anno[valid]] <- cb_chr$CYTOBAND[band_idx[valid]]
  }
  anno_df$CHR_CYTOBAND <- cytoband_vec

  # ---- DMR columns from bundled dmr_annotation ----
  dmr <- SEMseeker::dmr_annotation
  anno_df <- merge(anno_df, dmr, by = "PROBE", all.x = TRUE)

  # ---- Select final columns ----
  keep <- c(
    "PROBE", "CHR", "START", "END", tech,
    "GENE_BODY", "GENE_TSS200", "GENE_TSS1500", "GENE_1STEXON",
    "GENE_5UTR", "GENE_3UTR", "GENE_EXONBND", "GENE_WHOLE",
    "ISLAND_WHOLE", "ISLAND_N_SHORE", "ISLAND_S_SHORE",
    "ISLAND_N_SHELF", "ISLAND_S_SHELF",
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
