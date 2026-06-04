# NOTE: currently internal — could be re-exported once @examples are added
alphas_performance_pathway_analyser <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",
  significance = TRUE,disease_hpo, disease_description,keywords,stop_keywords,alphas,top=50,pathway_alpha=0.05, ...)
{
  #
  if(length(disease_hpo)>0)
    disease_original <- gsub("[:]","_",disease_hpo)
  else
    disease_original <- ""
  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  keys <- unique(ssEnv$keys_for_pathway)
  inference_details <- as.data.frame(inference_details)
  pvalue_column <- name_cleaning(pvalue_column)
  key_enrichment_format <- ssEnv$key_enrichment_format
  key_enrichment_format <- subset(key_enrichment_format, key_enrichment_format[,"label"]!="phenolyzer")

  columns_sorting <- c("SS_RANK","SS_RANK_ENRICHMENT","SS_RANK_FDR")
  for(column_sorting in columns_sorting)
  {
    for (i in 1:nrow(keys))
    {
      for (pt in 1:nrow(key_enrichment_format))
      {
        column_of_id <- key_enrichment_format[pt,"column_of_id"]
        column_of_adj_pvalue <- key_enrichment_format[pt,"column_of_adj_pvalue"]
        column_of_description <- key_enrichment_format[pt,"column_of_description"]
        column_of_enrichment <- key_enrichment_format[pt,"column_of_enrichment"]

        for (id in 1:nrow(inference_details))
        {
          inference_detail <- inference_details[id,]
          if(key_enrichment_format[pt,"type"]=="Pathway")
            path <- dir_check_and_create(ssEnv$result_folderPathway,c(key_enrichment_format[pt,"label"],name_cleaning(inference_detail$areas_sql_condition),name_cleaning(inference_detail$samples_sql_condition), name_cleaning(inference_detail$association_results_sql_condition)))
          else
            path <- dir_check_and_create(ssEnv$result_folderPhenotype,c(key_enrichment_format[pt,"label"],name_cleaning(inference_detail$areas_sql_condition),name_cleaning(inference_detail$samples_sql_condition), name_cleaning(inference_detail$association_results_sql_condition)))

          aggregated_patwhay_result_total <- data.frame()
          missed_keys <- data.frame()
          family_test <- inference_detail$family_test
          transformation_y <- as.character(inference_detail$transformation_y)
          independent_variable <- inference_detail$independent_variable
          covariates <- paste(inference_detail$covariates, collapse="_")

          file_prfx <- paste(independent_variable, transformation_y,family_test,covariates, pvalue_column,keys[i,"MARKER"], sep="_")


          for (alpha in alphas)
          {
            ssEnv$alpha <- alpha
            update_session_info(ssEnv)
            file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
              suffix="",
              pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
            file_name = file_path_build(path,file_name,"csv")
            # chech agin with disease without underscore
            if(!file.exists(file_name))
            {
              file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
                suffix="",
                pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
              file_name = file_path_build(path,file_name,"csv")
            }

            # check with without_signal at the end for pathfindR
            if(!file.exists(file_name))
            {
              file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
                suffix="",
                pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
              file_name = file_path_build(path,c(file_name,"without_signal"),"csv")
            }


            if (!file.exists(file_name))
            {
              # add a row if the file is missed it means non pathway was found !
              missed_keys <- plyr::rbind.fill(missed_keys, keys[i,])
              log_event("WARNING:", format(Sys.time(), "%a %b %d %X %Y") , " pathway_result ", file_name, " is missed !")
              next
            }

            pathway_result <- utils::read.csv2(file_name)

            if(nrow(pathway_result)==0)
              next

            # add fake FDR
            if(key_enrichment_format[pt,"label"]=="pathfindR")
              pathway_result$PVALUE_ADJ_ALL_FDR <-stats::p.adjust(pathway_result[,"highest_p"], method = "fdr")

            cols_to_check <- c(column_of_id,column_of_enrichment,column_of_description, column_of_adj_pvalue)
            # check column names contain the required columns
            if(!all(cols_to_check %in% colnames(pathway_result)))
            {
              log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y") , " column_of_id ", column_of_id, " is missed in ", file_name)
              next
            }

            pathway_result$MARKER <- keys[i,"MARKER"]
            pathway_result$FIGURE <- keys[i,"FIGURE"]
            # pathway_result$key <- paste(pathway_result$MARKER,pathway_result$FIGURE,pathway_result$AREA,pathway_result$SUBAREA,sep="_")


            pathway_result$by_keyword <- FALSE
            if(length(keywords)>0)
              pathway_result[,"by_keyword"] <- sapply(pathway_result[,column_of_description], function(x) any(grepl(paste(keywords, collapse = "|"), x, ignore.case = TRUE)))
            if(length(stop_keywords)>0)
              # remove item with the stop kewywords where by_keyword is TRUE
              pathway_result[pathway_result$by_keyword,"by_keyword"] <- sapply(pathway_result[pathway_result$by_keyword,column_of_description], function(x) !any(grepl(paste(stop_keywords, collapse = "|"), x, ignore.case = TRUE)))

            pathway_result$ALPHA <- ssEnv$alpha
            pathway_result <- enrichment_analysy_add_category(key_enrichment_format[pt,"label"],pathway_result)
            aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, pathway_result)
          }

          if(nrow(aggregated_patwhay_result_total)==0)
            next

          # TO DO: manage comparison of pathways by FIGURE
          aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, aggregated_patwhay_result_total$FIGURE =="HYPER_HYPO")

          # keep only first top taxonomies
          aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, aggregated_patwhay_result_total[,column_sorting] <= top)

          if(nrow(aggregated_patwhay_result_total) == 0)
            next

          # if(length(keywords)>0)
          #   aggregated_patwhay_result_total <- aggregated_patwhay_result_total[aggregated_patwhay_result_total[,column_of_adj_pvalue] < pathway_alpha,]

          if(nrow(aggregated_patwhay_result_total) == 0)
            next

          # if(any(aggregated_patwhay_result_total$by_keyword == TRUE))
          #   aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, by_keyword == TRUE)

          if(nrow(aggregated_patwhay_result_total) == 0)
            next

          path <- dir_check_and_create(ssEnv$result_folderChart,
            c("PATHWAYS",key_enrichment_format[pt,"label"],name_cleaning(inference_detail$areas_sql_condition),name_cleaning(inference_detail$samples_sql_condition),column_sorting))
          if (significance)
            marker_performance_pathway_plot(aggregated_patwhay_result_total, key_enrichment_format[pt,],file_prfx,path, "", performance_category="ALPHA", top=top, column_sorting)
        }
      }
    }
  }
}

