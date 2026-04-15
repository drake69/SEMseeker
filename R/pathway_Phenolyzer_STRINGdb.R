#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
pathway_Phenolyzer_STRINGdb <- function(study,
  statistic_parameter="",disease,
  adjust_per_area = FALSE, adjust_globally = FALSE,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
  inference_detail, significance = TRUE)
{

  #
  tmp <- tempdir()
  tempFolder <- paste(tmp,"/semseeker/",stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]"),sep="")
  pvalue_column <- name_cleaning(pvalue_column)

  # start_fresh <- FALSE
  # ssEnv <- init_env(result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
  ssEnv <- get_session_info()
  keys <- unique(ssEnv$keys_for_pathway)

  #check if optional package is installed
  if(!requireNamespace("STRINGdb", quietly = TRUE))
  {
    log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"),"STRINGdb package is not installed. Please install STRINGdb package to use this function.")
    return()
  }

  total_progress <- nrow(keys)*nrow(inference_detail)
  progress <- 0

  if(ssEnv$showprogress)
    progress_bar <- progressr::progressor(along = 1:total_progress)
  else
    progress_bar <- ""

  string_db <- STRINGdb::STRINGdb$new(version="11.5", species=9606, # 9606 is the taxonomy ID for Homo sapiens
    score_threshold=400, input_directory="")

  for(id in 1:nrow(inference_detail))
  {
    inference_detail <- inference_detail[id,]
    # foreach::foreach(i = 1:nrow(keys)) %dorng%
    for(i in 1:nrow(keys))
    {
      seq <- 0
      progress <- progress + 1

      if(exists("existing_db"))
        existing_keys <- unique(pathway_report[pathway_report$source==db,"key"])

      if(ssEnv$showprogress)
        progress_bar(sprintf("Searching for pathway using STRINGdb for %s",keys[i,]$COMBINED))

      key <- paste(keys[i,]$AREA,keys[i,]$SUBAREA,keys[i,]$MARKER,keys[i,]$FIGURE, sep="_")

      suffix <- ""
      if(statistic_parameter=="")
        suffix = "without_signal_"

      phenotype_analysis_name <- phenotype_analysis_name(inference_detail, keys[i,],prefix ="", suffix= suffix , pvalue_column, ssEnv$alpha, significance)
      path <- dir_check_and_create(ssEnv$result_folderPathway,c("Phenolyzer_STRINGdb",name_cleaning(inference_detail$areas_sql_condition), name_cleaning(inference_detail$samples_sql_condition), name_cleaning(inference_detail$association_results_sql_condition)))
      pathway_report_path <- file_path_build(path,phenotype_analysis_name,"csv")

      if(file.exists(pathway_report_path))
      {
        pp <- utils::read.csv2(pathway_report_path,stringsAsFactors = FALSE)
        if(nrow(pp)==0)
          next
        pathway_result_save(pp, pathway_report_path, "STRINGdb")
        next
      }
      #### START LOAD PHENOLYZER
      # load prioritized gene by phenolyzer
      base_path <- dir_check_and_create(ssEnv$result_folderPhenotype,c("phenolyzer",name_cleaning(inference_detail$areas_sql_condition)))
      phenotype_analysis_name <- phenotype_analysis_name( inference_detail = inference_detail,key = keys[i,], prefix="",suffix=paste("_", disease,"_report",sep=""), pvalue_column=pvalue_column, ssEnv$alpha, significance)
      path_phenolyzer <- dir_check_and_create(baseFolder = base_path, subFolders = "summary")
      phenotype_report_path <- file_path_build(path_phenolyzer,phenotype_analysis_name,"csv")

      if(!file.exists(phenotype_report_path))
      {
        log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " Phenotype report not found: ", phenotype_report_path)
        next
      }
      gene_set <- utils::read.csv2(phenotype_report_path,stringsAsFactors = FALSE)
      if(nrow(gene_set)==0)
        next
      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set: ",nrow(gene_set), " key: ", keys[i,])
      # log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set (ENTREZ): ",nrow(gene_set), " key: ", keys[i,])
      if (length(gene_set$Gene)==0)
        next

      #### END LOAD PHENOLYZER

      gene_set <- unique(gene_set)
      gene_set <- gene_set[order(gene_set[,"Rank"],decreasing = FALSE),]
      key <- paste(keys[i,]$AREA,keys[i,]$SUBAREA,keys[i,]$MARKER,keys[i,]$FIGURE, sep="_")

      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Number of genes in the gene set: ",nrow(gene_set), " key: ", keys[i,])

      if(nrow(gene_set)<2)
        next

      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Started STRINGdb analysis")

      gene_set <- gene_set[,c("Gene")]
      gene_set <- data.frame(pvalue=NA,logfc=NA,gene=gene_set,stringsAsFactors = FALSE)

      colnames(gene_set) <- c("pvalue","logfc","gene")
      # Map the gene list to STRING IDs
      mapped_genes <- string_db$map(gene_set, "gene")

      result_pathway <- string_db$get_enrichment(mapped_genes$STRING_id)
      if(nrow(result_pathway)==0)
        next
      result_pathway$key <- key

      # Assume total background number of genes (for Homo sapiens) is available
      # For example, let's assume the total background number of genes is 20000
      total_background_genes <- 20000

      # Calculate fold result_pathway
      result_pathway$expected_count <- (length(mapped_genes$STRING_id) * result_pathway$number_of_genes_in_background) / total_background_genes
      result_pathway$fold_enrichment <- result_pathway$number_of_genes / result_pathway$expected_count


      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Done STRINGdb analysis")


      if(exists("result_pathway"))
      {
        pathway_result_save(pp, pathway_report_path, "STRINGdb")
      }
    }

  }
}
