#' calculate scores for the same dataset at changing of quantile/bins
#'
#' @param inference_details data.frame where each row describes one inference run
#'   (columns: depth_analysis, family_test, transformation_y, covariates, independent_variable)
#' @param result_folder character string path to the root result folder
#' @param dest_folder character string path to the destination folder for aggregated output
#' @param folder_key character suffix appended to the study name to build per-q sub-folder names
#' @param pvalue_column name of the adjusted p-value column to use for significance ranking; default "PVALUE_ADJ_ALL_BH"
#' @param aggr_fun aggregation function name (as a string) used to summarise metrics across areas; default "mean"
#' @param significance logical; if TRUE, restrict comparison to significant results only
#' @param sql_conditions character vector of SQL condition strings used to subset samples
#' @param study character string identifying the study
#' @param qs integer vector of quantile/bin counts to evaluate; default c(2,4,8,16,32)
#' @param ... additional arguments passed to \code{init_env}
#'
#' @return A data.frame of performance metrics aggregated per marker and
#'   quantile/bin number; results are also written to CSV files on disk.
#'
#' @importFrom doRNG %dorng%
#' @importFrom doFuture %dofuture%
# NOTE: currently internal — could be re-exported once @examples are added
markers_performance_association_q_search <- function(inference_details, result_folder,dest_folder,folder_key, pvalue_column="PVALUE_ADJ_ALL_BH",
  aggr_fun ="mean",significance = TRUE,sql_conditions=c(),study, qs=c(2,4,8,16,32), ...)
{

  #
  # result_folder, maxResources = 90, parallel_strategy  = "multisession", ...
  #
  result_folder_temp <- paste0(result_folder,"/", study, folder_key, qs[1], sep="")
  ssEnv <- init_env( result_folder =  result_folder_temp, maxResources =  90, parallel_strategy  =  "sequential", start_fresh = FALSE)
  log_event("BANNER:", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will search evaluate the best quantile/bin post association analysis for project \n in ", ssEnv$result_folderData)

  markers <- unique(ssEnv$keys_areas_subareas_markers_figures$MARKER)
  markers <- sort(markers, decreasing = TRUE)

  model_metrics <- sort(c(ssEnv$model_metrics, pvalue_column))
  selected_figures <- unique(ssEnv$keys_markers_figures$FIGURE)
  inference_details <- as.data.frame(inference_details)


  # foreach::foreach(m = seq_along(markers)) %dorng%
  for (m in seq_along(markers))
  {
    marker <- markers[m]
    for (id in 1:nrow(inference_details))
    {
      #
      inference_detail <- inference_details[id,]
      depth_analysis <- inference_detail$depth_analysis
      family_test <- inference_detail$family_test
      transformation_y <- as.character(inference_detail$transformation_y)
      filtered_metrics <- metrics_filter(model_metrics, as.character(transformation_y))
      covariates <- paste0(inference_detail$covariates, collapse = "_")
      independent_variable <- inference_detail$independent_variable

      final <- data.frame()
      for (q in qs)
      {
        result_folder_temp <- dir_check_and_create(result_folder, paste0(study,folder_key,q))
        ssEnv <- init_env( result_folder =  result_folder_temp, maxResources =  90, parallel_strategy  =  "sequential", start_fresh = FALSE, verbosity=1)

        file_name <- inference_file_name(inference_detail, marker, ssEnv$result_folderInference,file_extension="csv", prefix = "", suffix="")
        # #
        if (!file.exists(file_name))
        {
          log_event("WARNING:", format(Sys.time(), "%a %b %d %X %Y") , "file ", file_name, " is missed !")
          next
        }

        # read the file
        file <- utils::read.csv2(file_name)


        # filter the metrics
        metrics_name_collect(file)

        file <- file[file$DEPTH == depth_analysis,]
        # file$COUNT <- 1
        file$KEY <- paste(file$FIGURE,file$AREA,file$SUBAREA, sep="_")
        # #
        file <- file[file$FIGURE %in% selected_figures,]

        file$Q <- q
        # bind the metrics
        final <- plyr::rbind.fill(final, file)
      }

      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Inference data loaded." )

      if (!exists("final"))
        next

      colnames(final) <- toupper(colnames(final))
      # filter using the sql condition
      final <- filter_sql(sql_conditions, final)

      if(nrow(final) == 0)
        next

      # count the number of rows for each marker, figure
      count <- aggregate(!(final$Q==""), by=list(final$Q), FUN=sum)
      # assign to each row the count of rows matching by marker and figure
      final$COUNT <- count$x[match(paste(final$Q),count$Group.1)]

      # remove entries with pvalue > alpha
      if (significance)
        final <- final[final[,pvalue_column] < as.numeric(ssEnv$alpha),]
      else
        final <- final[final[,pvalue_column] >= as.numeric(ssEnv$alpha),]
      # remove rows where p_value_column has NA
      final <- final[!is.na(final[,pvalue_column]),]
      if(nrow(final) == 0)
        next

      # count NA in the pvalue_column for each marker, figure
      count_na <- aggregate(final[,pvalue_column], by=list(final$Q), FUN=function(x) sum(is.na(x)))
      # assign to each row the count of NA in the pvalue_column matching by marker and figure
      final$COUNT_MISSED <- 0
      if(nrow(count_na) != 0)
      {
        count_na <- count_na[count_na$x >= 0,]
        # assign to each row the count of NA in the pvalue_column matching by marker and figure
        final$COUNT_MISSED <- count_na$x[match(paste(final$Q), count_na$Group.1)]
      }

      final$COUNT_MISSED <- (final$COUNT_MISSED / final$COUNT)

      # count number of rows for each marker and figure
      count <- aggregate(final$COUNT, by=list(final$Q), FUN=sum)
      # assign to each row the count of rows matching by marker and figure
      final$COUNT_SIGN <- count$x[match(paste(final$Q),count$Group.1)]

      final$COUNT_SIGN <- (final$COUNT_SIGN / final$COUNT)

      # check which exists in columns dataframe
      filtered_metrics <- unique(filtered_metrics[filtered_metrics %in% colnames(final)])
      if(depth_analysis > 2)
        filtered_metrics <- c(filtered_metrics, c("COUNT_SIGN","COUNT_MISSED"))
      filtered_metrics <- unique(sort(filtered_metrics))

      if(any("EFFECT_SIZE_MAGNITUDE" %in% filtered_metrics)){
        # replace large with 4 medium with 3 small with 2 and neglibible with 1
        final$EFFECT_SIZE_MAGNITUDE <- as.character(final$EFFECT_SIZE_MAGNITUDE)
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="large"] <- 4
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="medium"] <- 3
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="small"] <- 2
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="negligible"] <- 1
        final$EFFECT_SIZE_MAGNITUDE <- as.numeric(final$EFFECT_SIZE_MAGNITUDE)
      }

      filtered_metrics <- filtered_metrics[!is.na(filtered_metrics)]
      keys <- unique(final[,c("Q","FIGURE","AREA","SUBAREA")])
      keys <- unique(keys)

      area_subareas <- unique(keys[,c("AREA","SUBAREA")])
      if(ssEnv$showprogress)
        progress_bar <- progressr::progressor(along = 1:(nrow(area_subareas)*length(filtered_metrics)))
      else
        progress_bar <- ""

      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), "Starting score and plot creation." )

      for (as in 1:nrow(area_subareas))
      {
        subarea <- area_subareas[as,"SUBAREA"]
        area <- area_subareas[as, "AREA"]
        scores <- data.frame("Q"="","FIGURE"="","METRIC"="","SCORE"="")[-1,]
        figures <- unique(keys$FIGURE)
        # #
        for (f in seq_along(figures))
        {
          fig <- figures[f]
          for( m in seq_along(filtered_metrics)){
            metric <- filtered_metrics[m]
            # metric <- "MAE"
            final_temp <- final[final$FIGURE==fig,]
            final_temp <- final_temp[final_temp$AREA==area,]
            if(subarea=="ALL_SUBAREAS")
              final_temp <- final_temp[final_temp$SUBAREA!="WHOLE",]
            # group by Q and do mean of values of each group for metric
            final_temp <- final_temp[,colnames(final_temp) %in% (c("Q",metric))]
            # #
            if(nrow(final_temp) == 0)
              next


            final_temp <- metrics_ranking(metric = metric ,final_temp,column_to_rank =metric)
            scores_temp <- aggregate(final_temp[,"SCORE"], by = list(final_temp[,"Q"]), FUN = aggr_fun)
            scores_temp[,"FIGURE"] <- fig
            colnames(scores_temp) <- c("Q", "SCORE","FIGURE")
            scores <- rbind(scores, scores_temp)

            final_temp <- scores_temp[c("Q","SCORE")]
            colnames(final_temp) <- c("Q", "REBASED")

            if(ssEnv$showprogress)
              progress_bar(sprintf("Doing my job !"))

          }
        }

        #
        prfx <- paste(study, independent_variable, transformation_y,family_test, ssEnv$alpha ,covariates, area, subarea, marker, collapse="_")
        if(nrow(scores) != 0)
        {
          # score_max <- sum(scores$SCORE)
          # scores$SCORE <- round(100*scores$SCORE/score_max)
          # #
          # save scores
          path  <- dir_check_and_create(dest_folder,"MARKERS_PERFORMANCE")
          utils::write.csv2(scores, file = name_composer(path,"/",prfx,"_scores_", aggr_fun,".csv"), row.names = FALSE)

          scores$SCORE <- round(scores$SCORE,2)

          # aggregate scores by Q and sum RANK
          scores_agg <- aggregate(scores$SCORE, by = list(scores$Q), FUN = sum)
          # sort scores by SCORE descending
          scores_agg <- scores_agg[order(scores_agg$x, decreasing = TRUE),]
          colnames(scores_agg) <- c("Q","TOTAL")

          # create a pivot table with Q as rows and FIGURE as columns
          scores_agg_fig <- reshape2::dcast(scores, Q ~ FIGURE, value.var = "SCORE", fun.aggregate = sum)
          scores_agg_fig <- merge(scores_agg_fig,scores_agg, by="Q")
          # sort by SCORE descending
          scores_agg_fig <- scores_agg_fig[order(scores_agg_fig$TOTAL, decreasing = TRUE),]
          # save scores
          # utils::write.csv2(scores_agg_fig, file = name_composer(path,"/",prfx,"_scores_aggregated_",aggr_fun, ".csv"), row.names = FALSE)
        }


        #
        if(nrow(scores) != 0)
        {
          # score_max <- sum(scores$SCORE)
          # scores$SCORE <- round(100*scores$SCORE/score_max)
          # #
          # save scores
          path  <- dir_check_and_create(dest_folder,"MARKERS_PERFORMANCE")
          write.csv2(scores, file = name_composer(path,"/",prfx,"_scores_", aggr_fun,".csv"), row.names = FALSE)

          scores$SCORE <- round(scores$SCORE,2)

          # aggregate scores by Q and sum RANK
          scores_agg <- aggregate(scores$SCORE, by = list(scores$Q), FUN = sum)
          # sort scores by SCORE descending
          scores_agg <- scores_agg[order(scores_agg$x, decreasing = TRUE),]
          colnames(scores_agg) <- c("Q","TOTAL")

          # create a pivot table with Q as rows and FIGURE as columns
          scores_agg_fig <- reshape2::dcast(scores, Q ~ FIGURE, value.var = "SCORE", fun.aggregate = sum)
          scores_agg_fig <- merge(scores_agg_fig,scores_agg, by="Q")
          # sort by SCORE descending
          scores_agg_fig <- scores_agg_fig[order(scores_agg_fig$TOTAL, decreasing = TRUE),]
          # save scores
          utils::write.csv2(scores_agg_fig, file = name_composer(path,"/",prfx,"_scores_aggregated_",aggr_fun, ".csv"), row.names = FALSE)
        }
        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), "Areas, subareas done!" )

      }
      # }
    }
  }
  close_env()
}
