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
  group_column <- c("MARKER", "FIGURE", "AREA", "SUBAREA", "AREA_OF_TEST", "FAMILY_TEST",
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
