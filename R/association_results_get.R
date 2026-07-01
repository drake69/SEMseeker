# NOTE: naming convention — `adjustment_method` (singular) is used in internal
# helpers and result-retrieval functions where a single method string is expected.
# High-level public functions (e.g. enrichment_analysis()) use `adjustment_methods`
# (plural) because they accept a vector to iterate over multiple corrections.
association_results_get <- function (inference_detail, marker, adjust_per_area = FALSE, adjust_globally = FALSE,
  pvalue_column="PVALUE_ADJ_ALL_BH",adjustment_method = "BH", area ="GENE",
  omit_na = TRUE, significance = NULL)
{


  ssEnv <- core_get_session_info()
  resultFolder <- ssEnv$result_folderInference

  inferenceFile <- io_inference_file_name(inference_detail, marker, ssEnv$result_folderInference)
  if(adjust_per_area && adjust_globally)
  {
    core_log_event("ERROR: Can adjust per area or globbaly not both!", format(Sys.time(), "%a %b %d %X %Y"))
    stop()
  }

  if(!file.exists(inferenceFile))
  {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),  " Inference file does not exist: ", inferenceFile)
    return(data.frame())
  }


  results_inference <- utils::read.csv2(inferenceFile, row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
  colnames(results_inference) <- core_name_cleaning(colnames(results_inference))
  pvalue_column <- core_name_cleaning(pvalue_column)

  # AI-063: pvalue_columns in the user setup is a single vector applied to
  # every (inference_detail × marker) combination, but inference CSVs are
  # per-IV: a setup with IVs {TUMOUR_STAGE_N, BIOLOGICAL_RANK} writes
  # 'I_TUMOUR_STAGE_N_..._PVALUE_ADJ_ALL_FDR' to one CSV and
  # 'I_BIOLOGICAL_RANK_..._PVALUE_ADJ_ALL_FDR' to the other. Stopping the
  # whole run when one column happens to belong to the wrong IV kills the
  # entire enrichment pipeline on the first irrelevant lookup. Return an
  # empty result instead so the for-pvalue_columns loop in
  # enrichment_analysis() just skips this combination and keeps going.
  if((!pvalue_column %in% colnames(results_inference)))
  {
    core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " ", pvalue_column,
              " column does not exist in inference file (skipping): ",
              inferenceFile)
    return(data.frame())
  }

  # remove rows ehere pvalue_column is inf or -inf
  results_inference <- results_inference[!is.infinite(results_inference[,pvalue_column]),]

  sem_metrics_name_collect(results_inference)
  multiple_test_adj <- core_name_cleaning(ssEnv$multiple_test_adj)
  results_inference <- subset(results_inference,DEPTH==3)
  results_inference$SIGNIFICATIVE_ADJ <- apply(results_inference[, grepl(multiple_test_adj,colnames(results_inference))], 1, function(x) all(x < as.numeric(ssEnv$alpha)))
  results_inference$SIGNIFICATIVE <- apply(results_inference[, grepl("PVALUE", colnames(results_inference)) & !grepl(multiple_test_adj,colnames(results_inference))], 1, function(x) all(x < as.numeric(ssEnv$alpha)))
  # results_inference <- results_inference[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE","DEPTH")]

  results_inference[results_inference$AREA==area,"AREA_OF_TEST"] <- gsub(results_inference[results_inference$AREA==area,"AREA_OF_TEST"] , pattern="_", replacement="-")

  if(adjust_globally)
    results_inference[,pvalue_column] <-stats::p.adjust(results_inference[, pvalue_column], method = adjustment_method)

  adjustment_text <- "no_adjustment"
  if(adjust_per_area)
    adjustment_text <- "adjusted_per_area"
  if(adjust_globally)
    adjustment_text <- "adjusted_globally"
  # markers <- unique(keys$MARKER)

  areas <- unique(results_inference$AREA)
  if(adjust_per_area)
    for (a in seq_along(areas))
    {
      area <- areas[a]
      results_inference_area <- results_inference[results_inference$AREA==area,]
      results_inference_area[,pvalue_column] <-stats::p.adjust(results_inference_area[, pvalue_column], method = adjustment_method)
      results_inference[results_inference$AREA==area,] <- results_inference_area
    }

  results_inference <- subset(results_inference,AREA==area)

  # preserve only subareas selected
  results_inference <- results_inference[results_inference$SUBAREA %in% unique(ssEnv$keys_areas_subareas$SUBAREA),]
  if(!is.null(significance))
  {
    if (significance)
      results_inference <- subset(results_inference, results_inference[,pvalue_column] < as.numeric(ssEnv$alpha))
    else
      results_inference <- subset(results_inference, results_inference[,pvalue_column] >= as.numeric(ssEnv$alpha))
  }

  # remove where pvalue_column is NA
  results_inference <- results_inference[!is.na(results_inference[,pvalue_column]),]
  # remove where pvalue_column is -Inf or +Inf
  results_inference <- results_inference[results_inference[,pvalue_column] != -Inf,]
  results_inference <- results_inference[results_inference[,pvalue_column] != Inf,]

  # if(omit_na)
  #   results_inference <- na.omit(results_inference)

  results_inference <- filter_sql(inference_detail$areas_sql_condition, results_inference)

  if(nrow(results_inference) == 0)
    return(data.frame())

  entrez_ids <- rep(NA, length(results_inference$AREA_OF_TEST))
  tryCatch({
    entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = results_inference$AREA_OF_TEST, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  }, error = function(e) {

  })


  results_inference$ENTREZID <- entrez_ids
  #
  entrez_ids_not_na <- entrez_ids[!is.na(entrez_ids)]
  if ((length(entrez_ids_not_na) != length(results_inference$AREA_OF_TEST)))
  {
    lost_gene <- unique(results_inference$AREA_OF_TEST[is.na(entrez_ids)])
  }

  results_inference <- filter_sql(inference_details$association_results_sql_condition, results_inference)

  core_log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y")," inference file loaded:", inferenceFile)
  return(results_inference)
}
