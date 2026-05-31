# NOTE: currently internal — could be re-exported once @examples are added
markers_performance_pathway_analyser <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",
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

  for (a in alphas)
  {
    ssEnv$alpha <- a
    update_session_info(ssEnv)
    for (pt in 1:nrow(key_enrichment_format))
    {
      if(key_enrichment_format[pt,"label"]=="phenolyzer")
        disease <- disease_original
      else
        disease <- ""


      column_of_id <- key_enrichment_format[pt,"column_of_id"]
      column_of_adj_pvalue <- key_enrichment_format[pt,"column_of_adj_pvalue"]
      column_of_description <- key_enrichment_format[pt,"column_of_description"]
      column_of_enrichment <- key_enrichment_format[pt,"column_of_enrichment"]

      for (id in 1:nrow(inference_details))
      {
        # id <- 1
        # id <- 2
        inference_detail <- inference_details[id,]
        if(key_enrichment_format[pt,"type"]=="Pathway")
          path <- dir_check_and_create(ssEnv$result_folderPathway,c(key_enrichment_format[pt,"label"],name_cleaning(inference_detail$areas_sql_condition),name_cleaning(inference_detail$samples_sql_condition), name_cleaning(inference_detail$association_results_sql_condition)))
        else
          path <- dir_check_and_create(ssEnv$result_folderPhenotype,c(key_enrichment_format[pt,"label"],name_cleaning(inference_detail$areas_sql_condition),name_cleaning(inference_detail$samples_sql_condition), name_cleaning(inference_detail$association_results_sql_condition)))
        family_test <- inference_detail$family_test
        transformation_y <- as.character(inference_detail$transformation_y)
        independent_variable <- inference_detail$independent_variable
        covariates <- paste(inference_detail$covariates, collapse="_")

        file_prfx <- paste(independent_variable, transformation_y,family_test,covariates, pvalue_column,ssEnv$alpha, sep="_")

        aggregated_patwhay_result_total <- data.frame()
        missed_keys <- data.frame()
        for (i in 1:nrow(keys))
        {
          file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
            suffix=ifelse(disease=="","",paste("_",disease,sep="")),
            pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
          file_name = file_path_build(path,file_name,"csv")
          # chech agin with disease without underscore
          if(!file.exists(file_name))
          {
            file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
              suffix=ifelse(disease=="","",paste(disease,sep="")),
              pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
            file_name = file_path_build(path,file_name,"csv")
          }

          # check with without_signal at the end for pathfindR
          if(!file.exists(file_name))
          {
            file_name <- phenotype_analysis_name(inference_detail = inference_detail,key = keys[i,], prefix="",
              suffix=ifelse(disease=="","",paste(disease,sep="")),
              pvalue_column=pvalue_column, alpha = ssEnv$alpha, significance = significance)
            file_name = file_path_build(path,c(file_name,"without_signal"),"csv")
          }


          # #
          if (!file.exists(file_name))
          {
            # add a row if the file is missed it means non pathway was found !
            missed_keys <- plyr::rbind.fill(missed_keys, keys[i,])
            log_event("WARNING:", format(Sys.time(), "%a %b %d %X %Y") , " pathway_result ", file_name, " is missed !")
            next
          }

          if(key_enrichment_format[pt,"label"]=="phenolyzer")
          {
            # read the pathway_result
            pathway_result <- utils::read.csv2(file_name, dec=".")
          }
          else
            pathway_result <- utils::read.csv2(file_name)

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

          aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, pathway_result)
        }


        if(nrow(aggregated_patwhay_result_total)==0)
          next

        aggregated_patwhay_result_total <- enrichment_analysy_add_category(key_enrichment_format[pt,"label"],aggregated_patwhay_result_total)
        # TO DO: manage comparison of pathways by FIGURE
        aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, FIGURE =="HYPER_HYPO")



        # keep only first top taxonomies
        if(key_enrichment_format[pt,"label"]!="phenolyzer")
          aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, SS_RANK <= top)

        if(key_enrichment_format[pt,"label"]!="phenolyzer" & length(keywords)>0)
          aggregated_patwhay_result_total <- aggregated_patwhay_result_total[aggregated_patwhay_result_total[,column_of_adj_pvalue] <= pathway_alpha,]

        # remove rows where column_of_adj_pvalue is greater than alpha
        if(key_enrichment_format[pt,"label"]=="phenolyzer")
          aggregated_patwhay_result_total <- aggregated_patwhay_result_total[aggregated_patwhay_result_total[,column_of_adj_pvalue] > 0.33,]

        if(key_enrichment_format[pt,"label"]=="pathfindR")
          aggregated_patwhay_result_total <- aggregated_patwhay_result_total[aggregated_patwhay_result_total$support > 0.5,]

        # if(any(aggregated_patwhay_result_total$by_keyword == TRUE))
        #   aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, by_keyword == TRUE)

        if (!exists("aggregated_patwhay_result_total"))
          next
        if(nrow(aggregated_patwhay_result_total) == 0)
          next


        plot_path <- dir_check_and_create(ssEnv$result_folderChart,c("PATHWAYS",key_enrichment_format[pt,"label"],name_cleaning(inference_detail$areas_sql_condition),name_cleaning(inference_detail$samples_sql_condition)))
        if (significance & key_enrichment_format[pt,"label"]!="phenolyzer")
          marker_performance_pathway_plot(aggregated_patwhay_result_total, key_enrichment_format[pt,],file_prfx,plot_path, disease, performance_category="MARKER",top)

        # for each missed key add an empty row
        if(nrow(missed_keys)>0)
        {
          for (i in 1:nrow(missed_keys))
          {
            empty_row <- data.frame(column_of_id=NA, column_of_description =NA, column_of_adj_pvalue=NA, column_of_enrichment =NA, missed_keys[i,c("MARKER","FIGURE","AREA","SUBAREA")])
            empty_row$key <- paste(empty_row$MARKER,empty_row$FIGURE,empty_row$AREA,empty_row$SUBAREA,sep="_")
            colnames(empty_row) <- c(column_of_id, column_of_description, column_of_adj_pvalue, column_of_enrichment, "MARKER","FIGURE","AREA","SUBAREA","key")
            aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, empty_row)
          }
        }

        aggregated_patwhay_result <- aggregated_patwhay_result_total
        # save it
        utils::write.csv2(aggregated_patwhay_result, file_path_build(baseFolder =  path, detailsFilename =  paste(file_prfx,"_aggregated_patwhay_result", sep="_"), extension = "csv"))

        fdr <- aggregate(aggregated_patwhay_result[, column_of_adj_pvalue], by = list(aggregated_patwhay_result[,column_of_id]), FUN = mean)
        colnames(fdr) <- c( column_of_id, column_of_adj_pvalue)
        aggregated_patwhay_result <- unique(aggregated_patwhay_result[,c(column_of_id,"key",column_of_description)])
        aggregated_patwhay_result <- na.omit(aggregated_patwhay_result)
        categories <- unique(na.omit(aggregated_patwhay_result$key))
        # remove _
        categories <- gsub("_", " ", categories)
        if(length(categories)==1)
          next
        split <- split(aggregated_patwhay_result[,column_of_id], aggregated_patwhay_result$key)
        filename <- paste(path, "/",file_prfx,"_venn_diagramm_",key_enrichment_format[pt,"label"],".",ssEnv$plot_format,sep = "")


        # aggregated_patwhay_result$key <- gsub("_", " ", aggregated_patwhay_result$key)
        # change colname of column_of_id to column_of_id
        colnames(aggregated_patwhay_result)[which(colnames(aggregated_patwhay_result)==column_of_id)] <- "column_of_id_label"
        key_gene_set_pivot <- reshape2::dcast(aggregated_patwhay_result, column_of_id_label ~ key , value.var = "column_of_id_label", fun.aggregate = length)
        #change back colname of column_of_id to column_of_id
        colnames(aggregated_patwhay_result)[which(colnames(aggregated_patwhay_result)=="column_of_id_label")] <- column_of_id

        colnames(key_gene_set_pivot)[which(colnames(key_gene_set_pivot)=="column_of_id_label")] <- column_of_id

        # merge with geneSet name
        tt <- unique(aggregated_patwhay_result[,c(column_of_id,column_of_description)])
        tt <- merge(tt, fdr, by=column_of_id)
        key_gene_set_pivot <- merge(tt,key_gene_set_pivot, by=column_of_id)

        # add a column with the total number of gene sets
        key_gene_set_pivot$total <- rowSums(key_gene_set_pivot[,4:ncol(key_gene_set_pivot)])


        # foreach key add a column with column_of_enrichment
        kk <- unique(aggregated_patwhay_result$key)
        for (key in seq_along(kk))
        {
          pp <- aggregated_patwhay_result_total[aggregated_patwhay_result_total$key == kk[key],c(column_of_id,column_of_enrichment)]
          colnames(pp) <- c(column_of_id, paste("enrichment_",kk[key],sep = ""))
          key_gene_set_pivot <- merge(key_gene_set_pivot, pp, by=column_of_id, all.x = TRUE)
        }

        # add support to avaluate the quality of the gene set
        if(key_enrichment_format[pt,"label"]=="pathfindR")
          for (key in seq_along(kk))
          {
            pp <- aggregated_patwhay_result_total[aggregated_patwhay_result_total$key == kk[key],c(column_of_id,"support")]
            colnames(pp) <- c(column_of_id, paste("SUPPORT_",kk[key],sep = ""))
            key_gene_set_pivot <- merge(key_gene_set_pivot, pp, by=column_of_id, all.x = TRUE)
          }

        # key_gene_set_pivot <- key_gene_set_pivot[order(key_gene_set_pivot$total, decreasing = FALSE),]

        pivot_path <- dir_check_and_create(path, "marker_perfomance")
        filename <- paste(pivot_path, "/",file_prfx,"_pivot_",key_enrichment_format[pt,"label"],ifelse(disease=="","", paste("_", disease, sep="")), ".csv",sep = "")
        utils::write.csv2(key_gene_set_pivot, filename)

        # aggregate enrichment
        colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total)==column_of_id)] <- "column_of_id_label"
        enrichment <- reshape2::dcast(aggregated_patwhay_result_total, column_of_id_label ~ MARKER , value.var = column_of_enrichment, fun.aggregate = mean)
        colnames(enrichment)[which(colnames(enrichment)=="column_of_id_label")] <- column_of_id
        colnames(enrichment)[2:ncol(enrichment)] <- paste("ENRICHMENT_MEAN_",colnames(enrichment)[2:ncol(enrichment)],sep = "")
        #change back colname of column_of_id to column_of_id
        colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total)=="column_of_id_label")] <- column_of_id

        # aggregate fdr
        colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total)==column_of_id)] <- "column_of_id_label"
        fdr <- reshape2::dcast(aggregated_patwhay_result_total, column_of_id_label ~ MARKER , value.var = column_of_adj_pvalue, fun.aggregate = mean)
        colnames(fdr)[which(colnames(fdr)=="column_of_id_label")] <- column_of_id
        # add from 2 to end of columns FDR
        colnames(fdr)[2:ncol(fdr)] <- paste("FDR_MEAN_",colnames(fdr)[2:ncol(fdr)],sep = "")
        #change back colname of column_of_id to column_of_id
        colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total)=="column_of_id_label")] <- column_of_id

        colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total)==column_of_id)] <- "column_of_id_label"
        if(key_enrichment_format[pt,"label"]!="phenolyzer")
          # calculate min rank
          key_gene_set_pivot_summary <- reshape2::dcast(aggregated_patwhay_result_total, column_of_id_label ~ MARKER , value.var = "SS_RANK", fun.aggregate = min)
        else
          key_gene_set_pivot_summary <- reshape2::dcast(aggregated_patwhay_result_total, column_of_id_label ~ MARKER , value.var = "column_of_id_label", fun.aggregate = length)
        colnames(key_gene_set_pivot_summary)[which(colnames(key_gene_set_pivot_summary)=="column_of_id_label")] <- column_of_id
        colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total)=="column_of_id_label")] <- column_of_id

        # merge the total count
        key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary,key_gene_set_pivot[,c(column_of_id,"total")], by=column_of_id)
        key_gene_set_pivot_summary <- unique(merge(key_gene_set_pivot_summary,aggregated_patwhay_result[,c(column_of_id,column_of_description)], by=column_of_id, all.x = TRUE))

        # merge fdr and enrichment
        key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary, fdr, by=column_of_id, all.x=TRUE)
        key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary, enrichment, by=column_of_id, all.x=TRUE)

        # remove all columns with NA
        # key_gene_set_pivot_summary <- key_gene_set_pivot_summary[,colSums(is.na(key_gene_set_pivot_summary))<nrow(key_gene_set_pivot_summary)]

        # tryCatch(
        #   expr = {
        #
        #     if(ssEnv$openai_api_key !="")
        #     {
        #
        #       #
        #       # prompt_text <- paste(" I have a list of specific biological pathways and molecular functions, and I would like to understand if ",
        #       #   disease," is associated with each one.",
        #       #   " Can you identify which pathways or molecular functions from the list are linked to this disease and and provide the list of those linked ? ",
        #       #   " If none of the items on the list are linked, please return an empty string. \n", paste(batch, collapse = "\n"))
        #
        #
        #       # Function to create a prompt and get the response from OpenAI
        #       get_openai_response <- function(batch, disease) {
        #         #
        #         prompt_text <- paste("Give me back, just the list, of the passed items which are associated with ",
        #           disease,".","  If none are connected, give me back an empty string. \n", paste(batch, collapse = "\n"))
        #
        #         # prompt_text <- paste("I have a list of biological pathways or molecular function I'm passing to you. Give me back a list of those items which are associated with ",
        #         #   disease,".","  If none are connected, give me back an empty string. \n", paste(batch, collapse = "\n"))
        #
        #         # prompt_text <- paste(" I have a list of specific biological pathways and molecular functions, and I would like to understand if ",
        #         #   disease," is associated with each one.",
        #         #   " Can you identify which pathways or molecular functions from the list are linked to this disease and provide back the list of those linked ? ",
        #         #   " If none of the items on the list are linked, please return an empty string. \n", paste(batch, collapse = "\n"))
        #
        #         # prompt_text <- paste("I have a list of specific biological pathways and molecular functions, and I would like to understand if ",
        #         #   disease, " is associated with each one. ",
        #         #   "Can you identify which pathways or molecular functions from the list are linked to this disease and provide the list of those linked along with the likelihood of each association? ",
        #         #   "If none of the items on the list are linked, please return an empty string.\n",
        #         #   paste(batch, collapse = "\n"))
        #
        #         response <- openai::create_chat_completion(
        #           model = "gpt-4o",
        #           messages = list(
        #             list(role = "user", content = prompt_text)
        #           ),
        #           max_tokens = 1000,
        #           temperature = 0
        #         )
        #         response_text <- response$choices$message.content
        #         # split string by comma and remove whitespace
        #         response_vector <- unlist(strsplit(response_text, "\n"))
        #         rr <- batch %in% response_vector
        #         return(rr)
        #       }
        #
        #       # Process the items in batches
        #       batch_size <- 200  # Adjust batch size as needed
        #       descriptions <- key_gene_set_pivot_summary[, column_of_description]
        #       num_batches <- ceiling(length(descriptions) / batch_size)
        #       responses <- logical(length(descriptions))
        #
        #       for (i in 1:num_batches) {
        #         start_idx <- (i - 1) * batch_size + 1
        #         end_idx <- min(i * batch_size, length(descriptions))
        #         batch <- descriptions[start_idx:end_idx]
        #         responses[start_idx:end_idx] <- get_openai_response(batch, disease_description)
        #       }
        #
        #       #
        #       # # Create a prompt asking about the connection between each item and the phenotype
        #       # prompt <- paste("Give me back a vector of TRUE/FALSE to mark which of the passed items is connected with", disease, "\n", paste(key_gene_set_pivot_summary[,column_of_description], collapse = "\n"))
        #       # # Call the OpenAI API
        #       # response <- openai::create_completion(prompt, temperature = 0, max_tokens = 1)
        #       # # Extract the selected item from the response
        #       # selected_item <- response$choices[[1]]$text
        #       key_gene_set_pivot_summary[,"openai_suggested_item"] <- responses
        #     }
        #   },
        #   finally = {
        #   }
        # )
        #

        #

        # # calculate min excludinf Inf
        # key_gene_set_pivot_summary$min <- apply(key_gene_set_pivot_summary[,which(colnames(key_gene_set_pivot_summary) %in% unique(ssEnv$keys_markers_figures$MARKER))],
        #   1, function(x) min(x[!is.infinite(x)], na.rm = TRUE))
        #
        # # put in total the column name of the lowest value
        # key_gene_set_pivot_summary$total <- apply(key_gene_set_pivot_summary[,which(colnames(key_gene_set_pivot_summary) %in% unique(ssEnv$keys_markers_figures$MARKER))],
        #   1, function(x) colnames(key_gene_set_pivot_summary)[which(x==min(x, na.rm = TRUE))])

        key_gene_set_pivot_summary <- unique(merge(key_gene_set_pivot_summary, aggregated_patwhay_result_total[,c(column_of_id,"by_keyword")], by=column_of_id, all.x=TRUE))

        # put in total the column name of the lowest value excluding Inf
        key_gene_set_pivot_summary$total <- apply(key_gene_set_pivot_summary[,which(colnames(key_gene_set_pivot_summary) %in% unique(ssEnv$keys_markers_figures$MARKER))],
          1, function(x) {
            if (length(x[!is.infinite(x)]) == 1)
            {
              min_x <- x[!is.infinite(x)]
              sel_col <- which(x==min_x) + 1
            }
            else
            {
              min_x <- min(x[!is.infinite(x)], na.rm = TRUE)
              max_x <- max(x[!is.infinite(x)], na.rm = TRUE)
              if(min_x == max_x)
                return(NA)
              sel_col <- which(x==min_x) + 1
            }
            colnames(key_gene_set_pivot_summary)[sel_col]
          })




        if(key_enrichment_format[pt,"type"]=="Pathway")
          filename <- paste(ssEnv$result_folderPathway, "/",file_prfx,"_pivot_summary_",key_enrichment_format[pt,"label"], ifelse(disease=="","", paste("_", disease, sep="")) ,".csv",sep = "")
        else
          filename <- paste(ssEnv$result_folderPhenotype, "/",file_prfx,"_pivot_summary_",key_enrichment_format[pt,"label"], ifelse(disease=="","", paste("_", disease, sep="")) ,".csv",sep = "")
        utils::write.csv2(key_gene_set_pivot_summary, filename)



        # key_gene_set_pivot <- subset(key_gene_set_pivot, total < 3)[, colSums(key_gene_set_pivot[,4:ncol(key_gene_set_pivot)]) > 0,]
        # filename <- paste(path, "/",file_prfx,"_reduced_pivot_",key_enrichment_format[pt,"label"],".csv",sep = "")
        # write.csv2(key_gene_set_pivot, filename)

        # plot an heatmap
        # apcluster::heatmap(key_gene_set_pivot[,3:ncol(key_gene_set_pivot)], Rowv=NA, Colv=NA,
        #   col = colorRampPalette(c("white", "blue"))(100), scale="none", margins=c(5,10), main="Gene sets",
        #   xlab="Gene sets", ylab="Key", cexRow=0.5, cexCol=0.5, filename = paste(path, "/",file_prfx,"_heatmap_webgestalt.png",sep = ""))


        # #
        # Apply the function to the split data
        differences <- find_unique_gene_sets(split)
        # find differences
        # differences <- Reduce(setdiff, split)
        diff_df <- data.frame()
        #
        if (length(differences)>0)
        {
          filename <- paste(path, "/",file_prfx,"_differences_",key_enrichment_format[pt,"label"],".csv",sep = "")
          for (d in seq_along(differences))
          {
            tt2 <- merge(differences[d], names(differences[d]) )
            colnames(tt2) <- c(column_of_id,"KEY")
            diff_df <- plyr::rbind.fill(diff_df, tt2)
          }
          diff_df <- merge(diff_df, tt[,c(column_of_id, column_of_description)], by=column_of_id)
          utils::write.csv2(diff_df, filename)
        }
      }

      # remove all files ending with .log
      # files <- list.files(ssEnv$result_folderPathway, pattern = ".log")
      # for ( f in seq_along(files))
      # {
      #   file <- files[f]
      #   unlink(paste(ssEnv$result_folderPathway, "/", file, sep = ""), recursive = TRUE, force = TRUE)
      # }

    }
  }
}

