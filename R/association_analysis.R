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
#'     \item{depth_analysis}{Integer depth: \code{1} = sample level,
#'       \code{2} = type level (gene, DMR, CpG island),
#'       \code{3} = genomic area (TSS1550, WHOLE, TSS200, …).}
#'   }
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend; possible
#'   values: \code{"none"}, \code{"multisession"}, \code{"sequential"},
#'   \code{"multicore"}, \code{"cluster"} (default \code{"multicore"}).
#' @param start_fresh logical. If \code{TRUE}, delete previous inference
#'   results before running (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{init_env()}.
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

  ssEnv <- init_env(result_folder = result_folder, maxResources = maxResources,
    parallel_strategy = parallel_strategy, start_fresh = FALSE, ...)

  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
    " SEMseeker will perform the association analysys for project \n in ",
    ssEnv$result_folderData)

  if (start_fresh) unlink(ssEnv$result_folderInference, recursive = TRUE)
  dir_check_and_create(ssEnv$result_folderInference, c())

  localKeys <- ssEnv$keys_markers_figures

  deltaX_get()
  annotate_position_pivots()

  inference_details <- validate_inference_schema(unique(inference_details))

  for (z in seq_len(nrow(inference_details))) {
    start_time <- Sys.time()
    inference_detail <- inference_details[z, ]
    filter_p_value <- if (!is.null(inference_detail$filter_p_value))
      inference_detail$filter_p_value else TRUE

    log_inference_header(inference_detail)

    family_test <- split_and_clean(inference_detail$family_test)
    if (!validate_family_test(family_test)) next

    study_summary <- study_summary_get(inference_detail$samples_sql_condition)
    prep <- prepare_study_for_analysis(inference_detail, study_summary, family_test)
    if (is.null(prep)) next

    processed_items <- 0L
    last_results <- data.frame()
    last_filename <- NULL

    for (marker in unique(localKeys$MARKER)) {
      keys <- unique(localKeys[localKeys$MARKER == marker, ])
      fileNameResults <- inference_file_name(prep$inference_detail, marker,
        ssEnv$result_folderInference,
        prefix = ifelse(length(areas_selection) == 0, "",
          paste(areas_selection, "_", sep = "")))
      log_event("JOURNAL:", "Result saved into file:", fileNameResults, ".")

      # AI-040: skip sample-level (depth=1) and chr-level (depth=2) for
      # limma/voom families. They expect a per-area distribution to run
      # eBayes shrinkage on — depth=1 is a single-row fit (degenerate to
      # OLS) and depth=2 (TOTAL aggregate) mixes scales with depth=3 in
      # the same eBayes pool, contaminating the prior. Only depth=3
      # (per-probe / per-area) makes sense for these families.
      is_batch_family <- grepl("^(limma|voom)_", family_test)

      if (!is_batch_family) {
        d1 <- run_depth1_marker(prep, keys, family_test, fileNameResults,
          filter_p_value, ssEnv, ...)
        results <- d1$results
        processed_items <- processed_items + d1$processed_items
      } else {
        results <- data.frame()
        log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
                  " family_test='", family_test,
                  "': skipping DEPTH=1 (sample-level) — not meaningful for batch eBayes.")
      }

      if (prep$depth_analysis > 1) {
        dn <- run_depth_n_marker(prep, marker, family_test, fileNameResults,
          filter_p_value, ssEnv, selected_areas = areas_selection,
          results, start_time, processed_items, ...)
        results <- dn$results
        processed_items <- dn$processed_items
      }

      last_results  <- results
      last_filename <- fileNameResults

      # AI-061+ (2026-06-09): volcano plot for this marker right after the
      # CSV is finalised. One call per marker; volcano_plot_inference
      # splits internally by (AREA, SUBAREA) and writes one PNG per
      # combination under <result_folder>/Chart/VOLCANO/. Best-effort:
      # plot failure must not abort the analysis loop — log WARNING and
      # continue with the next marker.
      tryCatch(
        volcano_plot_inference(
          inference_detail = prep$inference_detail,
          result_folder    = ssEnv$result_folder,
          markers          = marker
        ),
        error = function(e) {
          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                    " volcano_plot_inference failed for marker '", marker,
                    "': ", conditionMessage(e))
        }
      )
    }

    finalize_job_results(last_results, prep$inference_detail, family_test,
      filter_p_value, last_filename, start_time, processed_items)
  }

  close_env()
}
