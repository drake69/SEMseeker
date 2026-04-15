pathway_result_save <- function(result_pathway, pathway_report_path, pathway_package)
{
  ssEnv <- get_session_info()
  if(nrow(result_pathway)!=0)
  {
    result_pathway <- enrichment_analysy_add_category(pathway_package,result_pathway)
    pvalue_column_adj <- ssEnv$key_enrichment_format[ssEnv$key_enrichment_format$label==pathway_package,"column_of_pvalue"]
    description_column <- ssEnv$key_enrichment_format[ssEnv$key_enrichment_format$label==pathway_package,"column_of_description"]
    # result_pathway <- result_pathway[, !(grepl("PVALUE_ADJ_", colnames(result_pathway))]
    col_p <- name_cleaning(paste0("PVALUE_ADJ_", ssEnv$multiple_test_adj))
    tryCatch({
      if(ssEnv$multiple_test_adj=="q")
        result_pathway[,col_p] <- qvalue::qvalue(result_pathway[,pvalue_column_adj], fdr.level = ssEnv$alpha, pi0.method="bootstrap", na.rm=TRUE)$qvalues
      else
        result_pathway[,col_p] <- stats::p.adjust(result_pathway[,pvalue_column_adj],method  =  ssEnv$multiple_test_adj)
      # sort by col_p
      result_pathway <- result_pathway[order(result_pathway[,col_p]),]
    }, error = function(e) {

    })
    result_pathway$PHENOTYPE <- grepl(study,result_pathway[,description_column], ignore.case = TRUE)
    # REMOVE COLUMNS with NAMES X, X.1 and X.2
    result_pathway <- result_pathway[,!grepl("^X$|^X\\.[0-9]+$", colnames(result_pathway))]
    utils::write.csv2(result_pathway, pathway_report_path, row.names = FALSE)
    rm(result_pathway)
  }
}
