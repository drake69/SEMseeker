assoc_analysis_save_results <- function(results=NULL,fileNameResults, family_test, filter_p_value, append=FALSE ){

  if(nrow(results)==0)
    return()

  ssEnv <- core_get_session_info()

  # C-06: stamp provenance columns so any downstream CSV stack retains build+tech
  genome_build_val <- if (!is.null(ssEnv$genome_build) && nzchar(ssEnv$genome_build))
    as.character(ssEnv$genome_build) else "hg19"
  tech_val <- if (!is.null(ssEnv$tech) && nzchar(ssEnv$tech))
    as.character(ssEnv$tech) else ""
  results$GENOME_BUILD <- genome_build_val
  results$TECH         <- tech_val

  utils::write.csv2(results,fileNameResults , row.names  =  FALSE)
  multiple_test_adj <- core_name_cleaning(ssEnv$multiple_test_adj)
  # there is a bug which mantain more family test in the same results file
  # so we need to filter the results
  #
  colnames(results) <- core_name_cleaning(colnames(results))
  results <- subset(results, FAMILY_TEST==as.character(family_test))

  # check if results is empty
  if(is.null(results))
    return()

  if(nrow(results)==0)
    return()

  if (!append)
  {
    results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
    results <- unique(results)

    pvalue_columns <- colnames(results)[grepl("PVALUE", colnames(results)) & !grepl("_ADJ", colnames(results))]

    # remove all existing column adjusted all pvalues
    results <- unique(results[,!grepl("_ADJ_ALL_", colnames(results))])

    # AI-257: PVALUE_ADJ is adjusted WITHIN THE KEY — the instances tested for
    # one measurement — and it is recomputed here because here is the only place
    # where every row of a key is together.
    #
    # assoc_apply_stat_model() adjusts what it has in hand, and what it has in
    # hand is one CHUNK: sem_run_depth_n_marker() splits a pivot at
    # ceiling(6e6 / ncol) rows, so a 485k-probe pivot on ~100 samples is nine
    # independent BH corrections instead of one. That made the family a memory
    # parameter rather than a statistical choice.
    #
    # The family is the key and not the whole file because the file now holds the
    # same area several times — once per aggregation asked for. Adjusting across
    # all of them would make the correction a function of how many ways you
    # looked rather than of how many hypotheses you tested: request MEAN, MEDIAN
    # and IQR of the same genes and every adjusted p roughly triples. The global
    # correction is still available in PVALUE_ADJ_ALL_<method>, and the two say
    # different things on purpose.
    results <- .assoc_adjust_within_key(results)

    if (exists("results") & length(pvalue_columns)>0)
    {
      for (p in seq_along(pvalue_columns))
      {
        col_p <- core_name_cleaning(paste0(pvalue_columns[p], "_ADJ_ALL_", multiple_test_adj))
        if(ssEnv$multiple_test_adj=="q")
          results[,col_p] <- qvalue::qvalue(results[,pvalue_columns[p]], fdr.level = ssEnv$alpha, pi0.method="bootstrap", na.rm=TRUE)$qvalues
        else
          results[,col_p] <- stats::p.adjust(results[,pvalue_columns[p]],method  =  ssEnv$multiple_test_adj)
        colnames(results) <- core_name_cleaning(colnames(results))
      }

      pvalue_adj_colname <- colnames(results)[grepl(multiple_test_adj,colnames(results))][1]

      if (nrow(results)>0)
        results <- results[order(results[,pvalue_adj_colname]),]

    }

    if(nrow(results)==0)
      return()

    results$DEPTH <- 3
    # replace NA of SUBAREA with TOTAL
    results[is.na(results$SUBAREA),"SUBAREA"] <- "TOTAL"
    results[results$SUBAREA=="SAMPLE","DEPTH"] <- 1
    selector <- grepl("TOTAL",results$AREA_OF_TEST)
    results[selector,"DEPTH"] <- 2
    # replace empty with NA
    results[results == ""] <- NA
    results[results == " "] <- NA
    # remove columns where all rows are NA
    results <- results[, colSums(is.na(results)) < nrow(results)]

    # check if exists at least a column with PVALUE
    if(!any(grepl("PVALUE", colnames(results))))
      return()
  }

  results$SIGNIFICATIVE_ADJ_ALL <- apply(as.data.frame(results[, grepl(multiple_test_adj,colnames(results))]), 1, function(x) all(x < as.numeric(ssEnv$alpha)))
  results$SIGNIFICATIVE <- apply(as.data.frame(results[, grepl("PVALUE", colnames(results)) & !grepl(multiple_test_adj,colnames(results))]), 1, function(x) all(x < as.numeric(ssEnv$alpha)))
  if(filter_p_value)
    results <- subset(results, SIGNIFICATIVE_ADJ)

  # remove duplicates based on MARKER   FIGURE  AREA    SUBAREA AREA_OF_TEST    FAMILY_TEST TRANSFORMATION_Y    PVALUE  R_MODEL
  # calculating the max of all others columns
  # C-06: include provenance columns in the grouping key so summarise() preserves them
  # AI-248: AGGREGATION is part of the identity. Without it the summarise(max)
  # below would fuse the median and the mean of the same scope into one row.
  group_column <- c("MARKER", "FIGURE", "AGGREGATION", "AREA", "SUBAREA", "AREA_OF_TEST", "FAMILY_TEST",
                    "TRANSFORMATION_Y", "R_MODEL", "TRANSFORMATION_X",
                    "INDEPENDENT_VARIABLE", "COVARIATES",
                    "GENOME_BUILD", "TECH")
  group_column <- group_column[group_column %in% colnames(results)]

  if(ncol(results[,!colnames(results) %in% group_column])>2)
  {
    # AI-061+ (2026-06-09): use base |> pipe (R 4.1+) instead of %>%.
    # The %>% reference was unresolved at runtime (no @importFrom magrittr)
    # and caused "could not find function %>%" mid-association on ewas v31.
    results <- results |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_column))) |>
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     ~ max(.x, na.rm = TRUE)),
                       .groups = 'drop')

    utils::write.csv2(results,fileNameResults , row.names  =  FALSE)
  }
}

#' Adjust p-values within the key, not within the chunk (internal)
#'
#' AI-257. The family over which an FDR is controlled has to be a statistical
#' choice. Two things had made it something else:
#'
#' \itemize{
#'   \item `assoc_apply_stat_model()` adjusts whatever `result_temp` it is
#'     handed, and at `SCOPE = INSTANCE` that is one chunk of the pivot — the
#'     family was the memory split;
#'   \item since the aggregation became an axis, one file holds the same area
#'     once per aggregation requested, so adjusting across the whole file makes
#'     the correction depend on how many ways the data were described.
#' }
#'
#' The family is therefore the identity key minus the instance: the areas tested
#' for one `(MARKER, FIGURE, SCOPE, AREA, SUBAREA, AGGREGATION)`. Asking for a
#' second aggregation no longer penalises the first, and the choice *between*
#' aggregations is a multiplicity to be handled by naming the aggregation in
#' advance — not by a correction, which cannot tell a planned comparison from an
#' opportunistic one.
#'
#' `PVALUE_ADJ_ALL_<method>` keeps the global view alongside it. The two answer
#' different questions and the analysis has to say which one it used.
#'
#' @param results the accumulated results of the run.
#' @return `results` with `PVALUE_ADJ` recomputed per family.
#' @keywords internal
#' @noRd
.assoc_adjust_within_key <- function(results) {

  if (is.null(results) || nrow(results) == 0 || !("PVALUE" %in% colnames(results)))
    return(results)

  key_cols <- intersect(c("MARKER", "FIGURE", "SCOPE", "AREA", "SUBAREA",
                          "AGGREGATION"),
                        colnames(results))
  if (length(key_cols) == 0)
    return(results)

  pvalues <- suppressWarnings(as.numeric(results$PVALUE))
  families <- do.call(paste, c(lapply(key_cols, function(cl)
    as.character(results[[cl]])), list(sep = "\r")))

  adjusted <- rep(NA_real_, length(pvalues))
  for (fam in unique(families)) {
    idx <- which(families == fam)
    adjusted[idx] <- stats::p.adjust(pvalues[idx], method = "BH")
  }

  results$PVALUE_ADJ <- adjusted
  results
}
