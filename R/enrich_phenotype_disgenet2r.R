# phenotype_disgenetplus2r <- function(study,
#   disease,score_limit,disgenetplus2r_user, disgenetplus2r_password,
#   pvalue = 0.05, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column="PVALUE_ADJ_ALL_BH",
#   inference_details,result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
# {
#
#   start_fresh <- FALSE
#   ssEnv <- core_init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = start_fresh, ...)
#
#   #check if optional package is installed
#   if(!requireNamespace("disgenetplus2r", quietly = TRUE))
#   {
#     core_log_event("disgenetplus2r package is not installed. Please install disgenetplus2r package to use this function with instruction: devtools::install_bitbucket('ibi_group/disgenetplus2r')")
#     return()
#   }
#
#   disgenet_api_key <- ""
#   tryCatch({
#     disgenet_api_key <- disgenetplus2r::get_disgenetplus_api_key(user = disgenetplus2r_user, password = disgenetplus2r_password,verbose = TRUE, warnings = TRUE)
#   }, error = function(e) {
#     core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " disgenetplus2r not working, internet problem or wrong credentials.")
#   })
#
#   if (is.null(disgenet_api_key))
#     return()
#
#   Sys.setenv(DISGENETPLUS_API_KEY= disgenet_api_key)
#
#
#   diseases <- enrich_disease_get(disease,"HPO")
#
#   if (is.null(nrow(diseases)))
#   {
#     core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " No disease found for ",disease)
#     return()
#   }
#
#   keys <- unique(ssEnv$keys_for_pathway)
#
#   if(ssEnv$showprogress)
#     progress_bar <- progressr::progressor(along = 1:(nrow(keys)))
#   else
#     progress_bar <- ""
#
#   seq <-0
#   for(i in 1:nrow(keys))
#   {
#     # i <-1
#     if(ssEnv$showprogress)
#       progress_bar(sprintf("Searching for disease using disgsenet2r: %s",keys[i,]$COMBINED))
#     key <- paste(keys[i,]$FIGURE,keys[i,]$MARKER,keys[i,]$AREA,keys[i,]$SUBAREA, sep="_")
#
#     enrich_phenotype_analysis_name <- enrich_phenotype_analysis_name( inference_detail = inference_details,key = keys[i,], prefix="",suffix=disease, pvalue_column=pvalue_column, pvalue, significance)
#     path <- io_dir_check_and_create(ssEnv$result_folderPhenotype,"disgenet")
#     phenotype_report_path <- io_file_path_build(path,enrich_phenotype_analysis_name,"csv")
#     if(file.exists(phenotype_report_path))
#       next
#
#     results_inference <- assoc_results_get(inference_detail,keys[i,"MARKER"],
#     pvalue, adjust_per_area, adjust_globally,pvalue_column,adjustment_method, significance = TRUE)
#
#     if (nrow(results_inference)==0)
#       next
#
#     if(keys[i,]$SUBAREA=="WHOLE")
#       gene_set <- results_inference[results_inference$SUBAREA==keys[i,]$SUBAREA,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
#     else
#       gene_set <- results_inference[results_inference$SUBAREA!=keys[i,]$SUBAREA,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
#
#     if(keys[i,]$FIGURE=="BOTH" | keys[i,]$FIGURE=="MEAN")
#       gene_set <- results_inference[results_inference$FIGURE==keys[i,]$FIGURE,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
#     else
#       gene_set <- results_inference[results_inference$FIGURE!="BOTH" & results_inference$FIGURE!="MEAN",c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column,"PVALUE"),]
#
#     if (nrow(results_inference)==0)
#       next
#
#
#     gene_set <- results_inference[,c("AREA_OF_TEST","STATISTIC_PARAMETER",pvalue_column)]
#     gene_set <- na.omit(gene_set)
#
#     if(nrow(gene_set)<2)
#       next
#
#     try(
#       {
#         entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.vector(gene_set$AREA_OF_TEST),column = "ENTREZID", keytype = "SYMBOL")
#         entrez <- unique(entrez)
#         entrez <- na.omit(entrez)
#
#         if (length(entrez)==0)
#           next
#
#         # res_enrich <- disgenetplus2r::gene2disease( gene =gene_set$AREA_OF_TEST, vocabulary = "ENTREZ", database = "ALL")
#         # gene_set <- as.vector(gene_set$AREA_OF_TEST)
#
#         res_enrich <- disgenetplus2r::gene2disease( gene =entrez, database = "ALL", score =c(0.5, 1), verbose  = TRUE)
#         output_temp <- disgenetplus2r::extract(res_enrich)
#
#         if(nrow(output_temp)==0)
#           next()
#
#         if(nrow(diseases)!=0)
#           output_temp$SEL <- output_temp$ID %in% diseases$DISEASE
#
#         seq <- seq + 1
#         output_temp$key <- key
#         output_temp$seq <- seq
#         output_temp$gene_count <- nrow(gene_set)
#         output_temp$order <- 1:nrow(output_temp)
#         output_temp$diseases <- paste(diseases$DISEASE, collapse = ";")
#         write.csv2(output_temp, phenotype_report_path)
#
#       }
#     )
#   }
# }
