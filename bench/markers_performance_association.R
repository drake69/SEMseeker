# NOTE: currently internal — could be re-exported once @examples are added
markers_performance_association <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",
  aggr_fun ="mean",significance = TRUE,sql_conditions=c(),alphas, ...)
{

  #
  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  unlink(path <- dir_check_and_create(ssEnv$result_folderChart, "MARKERS_PERFORMANCE"), recursive = TRUE)
  unlink( path  <- dir_check_and_create(ssEnv$result_folderInference,"MARKERS_PERFORMANCE"), recursive = TRUE)

  # for each marker read the file and extract the metrics and bind
  # them in a single dataframe
  markers <- unique(ssEnv$keys_areas_subareas_markers_figures$MARKER)
  # #

  #
  model_metrics <- sort(c(ssEnv$model_metrics, pvalue_column))
  selected_figures <- unique(ssEnv$keys_markers_figures$FIGURE)
  inference_details <- as.data.frame(inference_details)
  #

  for (a in alphas)
  {
    for (id in 1:nrow(inference_details))
    {

      ssEnv$alpha <- a
      update_session_info(ssEnv)
      inference_detail <- inference_details[id,]
      depth_analysis <- inference_detail$depth_analysis
      family_test <- inference_detail$family_test
      transformation_y <- as.character(inference_detail$transformation_y)
      if(transformation_y == "")
        transformation_y <- "none"
      filtered_metrics <- metrics_filter(model_metrics, as.character(transformation_y))
      covariates <- paste0(inference_detail$covariates, collapse = "_")
      independent_variable <- inference_detail$independent_variable

      final <- data.frame()
      for (m in seq_along(markers))
      {
        marker <- markers[m]
        # marker <- "DELTAQ"
        # message(marker)
        # inference_detail <- inference_details[1,]
        file_name <- inference_file_name(inference_detail, marker, ssEnv$result_folderInference,file_extension="csv", prefix = "", suffix="")
        # #
        if (!file.exists(file_name))
        {
          log_event("DEBUG: file ", file_name, " is missed !")
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

      # mantain only existng markers
      markers <- unique(final$MARKER)

      #
      # count the number of rows for each marker, figure
      count <- aggregate(!(final$MARKER==""), by=list(final$MARKER, final$FIGURE), FUN=sum)
      # assign to each row the count of rows matching by marker and figure
      final$COUNT <- count$x[match(paste(final$MARKER, final$FIGURE), paste(count$Group.1, count$Group.2))]


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
      count_na <- aggregate(final[,pvalue_column], by=list(final$MARKER, final$FIGURE), FUN=function(x) sum(is.na(x)))
      # assign to each row the count of NA in the pvalue_column matching by marker and figure
      final$COUNT_MISSED <- 0
      if(nrow(count_na) != 0)
      {
        count_na <- count_na[count_na$x >= 0,]
        # assign to each row the count of NA in the pvalue_column matching by marker and figure
        final$COUNT_MISSED <- count_na$x[match(paste(final$MARKER, final$FIGURE), paste(count_na$Group.1, count_na$Group.2))]
      }

      final$COUNT_MISSED <- (final$COUNT_MISSED / final$COUNT)

      # count number of rows for each marker and figure significative
      final_sign <- final[final[,pvalue_column] < as.numeric(ssEnv$alpha),]
      count <- aggregate(!(final_sign$MARKER==""), by=list(final_sign$MARKER, final_sign$FIGURE), FUN=sum)
      # assign to each row the count of rows matching by marker and figure
      final$COUNT_SIGN <- count$x[match(paste(final$MARKER, final$FIGURE), paste(count$Group.1, count$Group.2))]


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
      keys <- unique(final[,c("MARKER","FIGURE","AREA","SUBAREA")])
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
        # keys <- keys[keys$AREA==area & keys$SUBAREA==subarea,]
        # for (k in 1:nrow(keys))
        # {
        # create empty list to store plots
        plot_list <- list()
        plot_title_list <- list()
        row <- 0
        row_labels <- c()
        # key <- keys[k,]
        scores <- data.frame("MARKER"="","FIGURE"="","METRIC"="","SCORE"="")[-1,]
        figures <- unique(keys$FIGURE)
        # #
        for (f in seq_along(figures))
        {
          fig <- figures[f]
          row <- row  + 1
          col <- 0
          for( m in seq_along(filtered_metrics)){
            metric <- filtered_metrics[m]
            col <- col + 1
            # metric <- "MAE"
            final_temp <- final[final$FIGURE==fig,]
            final_temp <- final_temp[final_temp$AREA==area,]
            if(subarea=="ALL_SUBAREAS")
              final_temp <- final_temp[final_temp$SUBAREA!="WHOLE",]
            # group by MARKER and do mean of values of each group for metric
            final_temp <- final_temp[,colnames(final_temp) %in% (c("MARKER",metric))]
            # #
            if(nrow(final_temp) == 0)
              next
            ##############


            final_temp <- metrics_ranking(metric,final_temp, metric)
            scores_temp <- aggregate(final_temp[,"SCORE"], by = list(final_temp[,"MARKER"]), FUN = aggr_fun)
            scores_temp[,"FIGURE"] <- fig
            colnames(scores_temp) <- c("MARKER", "SCORE","FIGURE")
            scores <- rbind(scores, scores_temp)

            final_temp <- scores_temp[c("MARKER","SCORE")]
            colnames(final_temp) <- c("MARKER", "REBASED")
            # fill missed MARKER with metrics set to zero
            missed_markers <- setdiff(unique(final$MARKER), unique(final_temp$MARKER))
            if(length(missed_markers) > 0)
            {

              missed_markers <- data.frame("MARKER"=missed_markers, "REBASED"=rep(0,length(missed_markers)))
              final_temp <- rbind(final_temp, missed_markers)
            }

            ##############

            # final_temp <- aggregate(final_temp[,metric], by = list(final_temp[,"MARKER"]), FUN = ifelse(metric %in% c("COUNT_MISSED","COUNT_SIGN"),mean, aggr_fun))
            # colnames(final_temp) <- c("MARKER", "REBASED")
            # # #
            # # if any REBASE is not a number then skip
            # if(any(is.na(final_temp$REBASED)))
            #   next
            # scores <- metrics_ranking(metric,final_temp, scores, fig, "REBASED")

            bar_colors <- ssEnv$color_palette[seq_along(unique(final_temp$MARKER))]
            min_rebased <- min(final_temp$REBASED, na.rm = TRUE)
            max_rebased <- max(final_temp$REBASED, na.rm = TRUE)

            main_title <- ""
            if (row==1)
            {

              main_title <- gsub("_ADJ_ALL_","_ADJUSTED_",paste(metric))
              main_title <- gsub("_","\n",paste(main_title))
              main_title <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = main_title, size = 6, hjust = 0.5, vjust = 0.5) +
                ggplot2::theme_void()
              plot_title_list[[length(plot_title_list)+1]] <- main_title
            }
            y_title <- ""
            if (col==1)
              y_title <- fig

            # Create the plot
            ggp <- ggplot2::ggplot(final_temp, ggplot2::aes(x = MARKER, y = REBASED, fill = MARKER)) +
              ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
              ggplot2::scale_fill_manual(values = bar_colors) +  # Assign colors to bars
              ggplot2::labs(title = "", x = "", y = y_title) +
              ggplot2::theme_minimal() +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                legend.position = "none",
                plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::coord_cartesian(ylim = c(min_rebased, max_rebased))  # Set y-axis to start from min_rebased

            plot_list[[length(plot_list) + 1]] <- ggp

            if(ssEnv$showprogress)
              progress_bar(sprintf("Doing my job !"))

          }
        }

        # create a prefix for the file name
        prfx <- paste( independent_variable, transformation_y,family_test, ssEnv$alpha ,covariates, area, subarea, paste0(markers,collapse="_"))
        if(length(plot_list) != 0)
        {

          ggt <- gridExtra::grid.arrange(grobs = plot_title_list, nrow = 1)
          path <- dir_check_and_create(ssEnv$result_folderChart, "MARKERS_PERFORMANCE")
          file_name <- inference_file_name(inference_detail, prfx,path,file_extension=ssEnv$plot_format, prefix = "_TITLE_", suffix=aggr_fun)
          # save the panel
          ggplot2::ggsave(file = file.path(file_name), ggt, width = 4960, height = 400, units = "px")

          # build  a panel from plot list
          gge <- gridExtra::grid.arrange(grobs = lapply(plot_list, ggplot2::ggplotGrob), ncol = length(filtered_metrics))
          path <- dir_check_and_create(ssEnv$result_folderChart, "MARKERS_PERFORMANCE")
          file_name <- inference_file_name(inference_detail, prfx,path,file_extension=ssEnv$plot_format, prefix = "", suffix=aggr_fun)
          # save the panel
          ggplot2::ggsave(file = file.path(file_name), gge, width = 4960, height = 2480, units = "px")

        }

        if(nrow(scores) != 0)
        {
          # score_max <- sum(scores$SCORE)
          # scores$SCORE <- round(100*scores$SCORE/score_max)
          # #
          # save scores
          path  <- dir_check_and_create(ssEnv$result_folderInference,"MARKERS_PERFORMANCE")
          utils::write.csv2(scores, file = paste0(path,"/",prfx,"_scores_", aggr_fun,".csv"), row.names = FALSE)


          # aggregate scores by MARKER and sum RANK
          scores_agg <- aggregate(scores$SCORE, by = list(scores$MARKER), FUN = sum)
          # sort scores by SCORE descending
          scores_agg <- scores_agg[order(scores_agg$x, decreasing = TRUE),]
          colnames(scores_agg) <- c("MARKER","TOTAL")

          # create a pivot table with MARKER as rows and FIGURE as columns
          scores_agg_fig <- reshape2::dcast(scores, MARKER ~ FIGURE, value.var = "SCORE", fun.aggregate = sum)
          scores_agg_fig <- merge(scores_agg_fig,scores_agg, by="MARKER")
          # sort by SCORE descending
          scores_agg_fig <- scores_agg_fig[order(scores_agg_fig$TOTAL, decreasing = TRUE),]
          # save scores
          fname <- paste0(path,"/",prfx,"_scores_aggregated_",aggr_fun, ".csv")
          save_latex_table(scores_agg_fig, fname, "Scores post association-analysis per each marker.")
          utils::write.csv2(scores_agg_fig, file = fname, row.names = FALSE)
        }

        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), "Areas, subareas done!" )

      }
      # }
    }
  }
  close_env()
}
