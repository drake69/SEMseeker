#' Association analysis of SEMseeker results
#'
#' Run statistical association models between SEM metrics and a phenotype
#' variable. Supports group tests (Wilcoxon, t-test), GLM families (gaussian,
#' poisson, binomial), quantile regression, correlations (Pearson, Kendall,
#' Spearman), and multi-covariate formulas (e.g.
#' \code{MUTATIONS_* ~ covariate1 + covariate2}).
#'
#' @param inference_details data.frame. Each row defines one analysis run.
#'   Required columns:
#'   \describe{
#'     \item{independent_variable}{Sample sheet column used as grouping /
#'       covariate variable.}
#'     \item{family_test}{Statistical model: \code{"wilcoxon"},
#'       \code{"stats::t.test"}, \code{"gaussian"}, \code{"poisson"},
#'       \code{"binomial"}, \code{"pearson"}, \code{"kendall"},
#'       \code{"spearman"}, or quantile regression as
#'       \code{"quantreg_<tau>_<runs>"} (e.g. \code{"quantreg_0.25_2000"}).}
#'     \item{transformation_y}{Transformation applied to the dependent variable:
#'       \code{"none"}, \code{"scale"}, \code{"log"}, \code{"log2"},
#'       \code{"log10"}, \code{"exp"}, or
#'       \code{"quantile_<n>"} (e.g. \code{"quantile_3"}).}
#'     \item{marker}{SEM metric column prefix (e.g. \code{"DELTARP"},
#'       \code{"MUTATIONS"}).}
#'     \item{depth_analysis}{Whether the region class is collapsed or kept open.
#'       \code{1} tests one number per sample (scope \code{SAMPLE}); anything
#'       above tests one row per instance of the class — a gene, an island, a
#'       cytoband, a probe (scope \code{INSTANCE}). It used to be documented as a
#'       three-step ladder (sample, type, genomic area), which described a
#'       containment that the code never implemented and that the region classes
#'       do not have: a TSS200 window and an open-sea stretch refine neither each
#'       other nor anything between them, so they cannot be placed on one line.
#'       What the artefact covers is said by \code{AREA} and \code{SUBAREA}; this
#'       value survives as a label on the result rows.}
#'     \item{aggregation}{Required. How the positions are reduced to the one
#'       number the model is fitted on: \code{"SUM"}, \code{"MEAN"},
#'       \code{"MEDIAN"}, \code{"VARIANCE"}, \code{"IQR"}, \code{"MODELOW"} or
#'       \code{"MODEHIGH"}. While every marker admitted exactly one operator
#'       this could stay implicit; a scope now carries several, so the request
#'       has to name the one it wants. Which are admissible depends on the
#'       marker: a count carries \code{SUM} (the burden) and \code{MEAN} (the
#'       density, the form comparable across regions of different size), while
#'       its median and IQR are degenerate; the two modes exist only for the
#'       signal on the beta scale. A request no marker of the run admits is
#'       dropped with a warning naming it, not answered with an empty result.}
#'     \item{scopes}{Optional. Which region classes to collapse to one number
#'       per sample, several separated by \code{"+"} (default \code{"SAMPLE"},
#'       meaning no restriction). A class such as \code{"GENE_TSS1500"} tests the
#'       burden restricted to the positions of that class, and the result rows
#'       carry \code{SCOPE = "SAMPLE"} with \code{AREA = "GENE"} and
#'       \code{SUBAREA = "TSS1500"} — the coordinates of the taxonomy, not the
#'       scope name squashed into one of them. The per-instance artefacts of the
#'       run are tested regardless, and carry \code{SCOPE = "INSTANCE"}.
#'
#'       The artefact is built on the way in if it does not exist yet (AI-255),
#'       so a class no previous run foresaw costs one scan of the position pivot
#'       rather than a rerun — provided this call registers it in \code{areas}
#'       and \code{subareas}.}
#'   }
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend; possible
#'   values: \code{"none"}, \code{"multisession"}, \code{"sequential"},
#'   \code{"multicore"}, \code{"cluster"} (default \code{"multicore"}).
#' @param start_fresh logical. If \code{TRUE}, delete previous inference
#'   results before running (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{core_init_env()}.
#'
#' @return Invisibly \code{NULL}. Inference result CSV files are written to
#'   the \code{Inference/} sub-folder of \code{result_folder}, one file per
#'   marker/area/family combination defined in \code{inference_details}.
#' @importFrom doRNG %dorng%
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' association_analysis(
#'   inference_details = data.frame(
#'     independent_variable = "Sample_Group",
#'     family_test          = "wilcoxon",
#'     transformation_y     = "none",
#'     marker               = "DELTARP",
#'     areas                = "GENE"
#'   ),
#'   result_folder     = "~/semseeker_results/",
#'   multiple_test_adj = "BH"
#' )
#' }
#' @export
association_analysis <- function(inference_details, result_folder, maxResources = 90,
  parallel_strategy = "multicore", start_fresh = FALSE, ...) {

  arguments <- list(...)
  areas_selection <- c()
  if (!is.null(arguments[["areas_selection"]])) {
    areas_selection <- arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- core_init_env(result_folder = result_folder, maxResources = maxResources,
    parallel_strategy = parallel_strategy, start_fresh = FALSE, ...)

  core_log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
    " SEMseeker will perform the association analysys for project \n in ",
    ssEnv$result_folderData)

  if (start_fresh) unlink(ssEnv$result_folderInference, recursive = TRUE)
  io_dir_check_and_create(ssEnv$result_folderInference, c())

  localKeys <- ssEnv$keys_markers_figures

  sem_deltaX_get()
  anno_annotate_position_pivots()

  inference_details <- assoc_validate_inference_schema(unique(inference_details))
  # AI-248: shape first, then meaning. Refuse an aggregation that no marker of
  # this run admits before any result is written — checking it inside the
  # per-marker loop would surface the mistake after part of the output exists.
  inference_details <- assoc_validate_aggregation(inference_details)

  for (z in seq_len(nrow(inference_details))) {
    start_time <- Sys.time()
    inference_detail <- inference_details[z, ]
    filter_p_value <- if (!is.null(inference_detail$filter_p_value))
      inference_detail$filter_p_value else TRUE

    core_log_inference_header(inference_detail)

    family_test <- util_split_and_clean(inference_detail$family_test)
    if (!assoc_validate_family_test(family_test)) next

    # AI-255: the models read artefacts, not columns — assoc_run_marker() opens
    # the pivot for every key, collapsed or not. So what this needs from the
    # sample sheet is the phenotype and the covariates, and joining the
    # per-sample statistics onto it would build artefacts nobody then reads:
    # io_feature_colname() has exactly one caller left, the composer inside
    # sem_study_summary_get(), and nothing reads those names back.
    #
    # The join is still done when the request names a feature the plain sheet
    # does not have — adjusting for the global burden is a legitimate thing to
    # ask — but it is no longer paid for on every run by default.
    study_summary <- sem_study_summary_get(inference_detail$samples_sql_condition,
                                           with_sample_stats = FALSE)
    wanted_cols <- c(gsub(" ", "", as.character(inference_detail$independent_variable)),
                     util_split_and_clean(inference_detail$covariates))
    wanted_cols <- wanted_cols[nzchar(wanted_cols) & !is.na(wanted_cols)]
    if (!is.null(study_summary) && !all(wanted_cols %in% colnames(study_summary))) {
      core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                " The request names ", paste(setdiff(wanted_cols, colnames(study_summary)),
                                             collapse = ", "),
                ", which the sample sheet does not carry: joining the per-sample ",
                "features as well.")
      study_summary <- sem_study_summary_get(
        inference_detail$samples_sql_condition,
        regions = util_split_and_clean(inference_detail$scopes))
    }
    prep <- sem_prepare_study_for_analysis(inference_detail, study_summary, family_test)
    if (is.null(prep)) next

    processed_items <- 0L
    last_results <- data.frame()
    last_filename <- NULL

    for (marker in unique(localKeys$MARKER)) {
      keys <- unique(localKeys[localKeys$MARKER == marker, ])
      fileNameResults <- io_inference_file_name(prep$inference_detail, marker,
        ssEnv$result_folderInference,
        prefix = ifelse(length(areas_selection) == 0, "",
          paste(areas_selection, "_", sep = "")))
      core_log_event("JOURNAL:", "Result saved into file:", fileNameResults, ".")

      # AI-255: one road. There used to be two calls here, chosen by
      # depth_analysis, because the collapsed artefact and the per-instance one
      # had different shapes — a table of columns against a pivot of rows. They
      # have the same shape now, so a model handed a row does not know, and has
      # no reason to ask, whether the key of that row is a gene symbol or
      # PROBE_WHOLE. It fits. The scope travels in the key; the batch-family
      # exclusion travels with it (see .assoc_marker_keys).
      dn <- assoc_run_marker(prep, marker, family_test, fileNameResults,
        filter_p_value, ssEnv, selected_areas = areas_selection,
        data.frame(), start_time, processed_items, ...)
      results <- dn$results
      processed_items <- dn$processed_items

      last_results  <- results
      last_filename <- fileNameResults

      # AI-061+ (2026-06-09): volcano plot for this marker right after the
      # CSV is finalised. One call per marker; assoc_volcano_plot_inference
      # splits internally by (AREA, SUBAREA) and writes one PNG per
      # combination under <result_folder>/Chart/VOLCANO/. Best-effort:
      # plot failure must not abort the analysis loop — log WARNING and
      # continue with the next marker.
      tryCatch(
        assoc_volcano_plot_inference(
          inference_detail = prep$inference_detail,
          result_folder    = ssEnv$result_folder,
          markers          = marker
        ),
        error = function(e) {
          core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                    " assoc_volcano_plot_inference failed for marker '", marker,
                    "': ", conditionMessage(e))
        }
      )
    }

    util_finalize_job_results(last_results, prep$inference_detail, family_test,
      filter_p_value, last_filename, start_time, processed_items)
  }

  core_close_env()
}
