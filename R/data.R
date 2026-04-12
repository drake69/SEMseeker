#' dmr_annotation
#'
#' Differentially methylated region (DMR) annotations for CpG probes on
#' Illumina methylation arrays.  Contains only probes that overlap at least one
#' known imprinted or disease-associated DMR.  This lightweight table
#' (~1,600 rows) supplements the Bioconductor array annotation packages used by
#' \code{\link{probe_annotation_build}}, which do not carry DMR-level
#' annotations.
#'
#' @format A data frame with three columns:
#' \describe{
#'   \item{PROBE}{Illumina probe identifier (e.g. \code{"cg00000924"}).}
#'   \item{DMR_WHOLE}{Genomic region label of the enclosing DMR
#'     (e.g. \code{"KCNQ1OT1:TSS-DMR"}).}
#'   \item{DMR_DMR}{Fine-grained DMR label, identical to \code{DMR_WHOLE} for
#'     most probes.}
#' }
"dmr_annotation"


#' cytoband_hg19
#'
#' Cytogenetic band coordinates for the hg19 human genome assembly.
#' Used by \code{\link{probe_annotation_build}} to assign a \code{CHR_CYTOBAND}
#' label (e.g. \code{"q12.2"}) to each CpG probe based on its chromosomal
#' position.  Band boundaries are derived from probe positions in the full
#' Illumina EPIC annotation and cover all autosomes and sex chromosomes.
#'
#' @format A data frame with 829 rows and four columns:
#' \describe{
#'   \item{CHR}{Chromosome identifier without \code{"chr"} prefix
#'     (e.g. \code{"1"}, \code{"X"}).}
#'   \item{START}{Approximate start position of the cytogenetic band (bp).}
#'   \item{END}{Approximate end position of the cytogenetic band (bp).}
#'   \item{CYTOBAND}{ISCN band label (e.g. \code{"q12.2"}, \code{"p36.33"}).}
#' }
"cytoband_hg19"


#' metrics_properties
#'
#' Metadata table describing the statistical properties of each SEM metric.
#' Used internally to determine ranking direction and scaling behaviour during
#' association analysis and pathway enrichment.
#'
#' @format A data frame with one row per metric and columns including
#'   \code{Metric} (metric name), \code{Higher_the_Better} (logical: whether
#'   higher values indicate stronger signal), and \code{Affected_by_Scaling}
#'   (logical: whether the metric is affected by data transformation).
"metrics_properties"


#' ssEnv
#'
#' Internal session environment object persisted between SEMseeker analysis
#' steps. Stores runtime parameters (result folder paths, technology flag,
#' alpha threshold, etc.) set by \code{\link{init_env}} and retrieved by
#' \code{get_session_info()}.
#'
#' @format An \code{environment} containing named slots for session-level
#'   analysis parameters. Users should not modify this object directly;
#'   use \code{\link{init_env}} and \code{set_env_variable()} instead.
"ssEnv"
