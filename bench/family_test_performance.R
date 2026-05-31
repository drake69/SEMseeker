family_test_performance <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH", ...)
{
  inference_details <- as.data.frame(inference_details)
  # all families must have the same scale
  transformations <- unique(inference_details$transformation_y)
  depths <- unique(inference_details$depth_analysis)
  # ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE)
  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  # for each family test, read the file and extract the metrics and bind
  # them in a single dataframe
  # 
  markers <- unique(ssEnv$keys_markers_figures$MARKER)
  model_metrics <- sort(c(ssEnv$model_metrics, pvalue_column))
  selected_figures <- unique(ssEnv$keys_markers_figures$FIGURE)

  for ( d in seq_along(depths))
  {
    depth <- depths[d]
    for (i in seq_along(transformations))
    {
      transformation_y <-transformations[i]
      filtered_metrics <- metrics_filter(model_metrics, transformation_y)

      inference_details_selection <- inference_details[inference_details$transformation_y==transformation_y & inference_details$depth_analysis==depth,]
      # check different families.test
      families.test <- unique(inference_details_selection$family_test)

      if (length(families.test)!=nrow(inference_details_selection))
        stop("All families.test must be different.")


      for (m in seq_along(markers))
      {
        marker <- markers[m]
        # marker <- "DELTAQ"
        for (i in 1:nrow(inference_details_selection))
        {
          # i <- 1
          inference_detail <- inference_details_selection[i,]
          family_test <- inference_details_selection[i,"family_test"]
          file_name <- inference_file_name(inference_detail, marker, ssEnv$result_folderInference,file_extension="csv", prefix = "", suffix="")
          if(!file.exists(file_name))
          {
            log_event("DEBUG: File does not exist: ", file_name, " skipping")
            next
          }
          # read the file
          file <- utils::read.csv2(file_name)
          file <- file[file$DEPTH == depths,]
          # file$FAMILY_TEST <- paste(family_test, inference_detail$independent_variable, sep="_")
          metrics_name_collect(file)
          file$KEY <- paste(file$FIGURE,file$AREA,file$SUBAREA, sep="_")
          file <- file[file$FIGURE %in% selected_figures,]
          # bind the metrics
          if (!exists("final"))
            final <- file
          else
            final <- plyr::rbind.fill(final, file)
        }
      }

      # remove entries with pvalue > alpha
      final <- final[final[,pvalue_column] < as.numeric(ssEnv$alpha),]

      # remove rows where p_value_colimn has NA
      final <- final[!is.na(final[,pvalue_column]),]

      # remove columns where any NA
      final <- final[,colSums(is.na(final))==0]



      # check which exists in columns dataframe
      filtered_metrics <- unique(filtered_metrics[filtered_metrics %in% colnames(final)])

      if(any("EFFECT_SIZE_MAGNITUDE" %in% filtered_metrics)){
        # replace large with 4 medium with 3 small with 2 and neglibible with 1
        final$EFFECT_SIZE_MAGNITUDE <- as.character(final$EFFECT_SIZE_MAGNITUDE)
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="large"] <- 4
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="medium"] <- 3
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="small"] <- 2
        final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="negligible"] <- 1
        final$EFFECT_SIZE_MAGNITUDE <- as.numeric(final$EFFECT_SIZE_MAGNITUDE)
      }

      for (m in seq_along(markers))
      {
        marker <- markers[m]
        if (nrow(final[final$MARKER==marker,]) == 0)
        {
          log_event("DEBUG: No data for marker: ", marker, " skipping")
          next
        }
        # marker <- "DELTAQ"
        # creat empty list to store plots
        plot_list <- list()
        # loop over metrics
        for ( m in seq_along(filtered_metrics))
        {
          filtered_metric <- filtered_metrics[m]
          # filtered_metric <- "MAE"
          # create a plot
          p <- ggplot2::ggplot(final[final$MARKER==marker,c("FAMILY_TEST",filtered_metric)], ggplot2::aes_string(x="FAMILY_TEST", y=filtered_metric, fill="FAMILY_TEST")) +
            ggplot2::geom_boxplot() +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
              legend.position = "none") +
            ggplot2::labs(title = "") +
            ggplot2::ylab(filtered_metric) +
            ggplot2::xlab("") +
            ggplot2::scale_fill_manual(values = ssEnv$color_palette)
          # store the plot in the list
          plot_list[[filtered_metric]] <- p
        }

        # build  a panel from plot list
        gge <- gridExtra::grid.arrange(grobs = lapply(plot_list, ggplot2::ggplotGrob), ncol = length(filtered_metrics))

        path <- dir_check_and_create(ssEnv$result_folderChart, "FAMILY.TEST_PERFORMANCE")
        # file_name <- inference_file_name(inference_detail, marker,path ,file_extension=ssEnv$plot_format, prefix = "", suffix="")
        file_name <- file.path(path, paste0(transformation_y,"_", marker,"_",paste0(families.test, collapse="_"),"_FAMILY_TEST_PERFORMANCE.png"))
        # save the panel
        ggplot2::ggsave(file = file.path(file_name), gge, width = 2048, height = 1024, units = "px")
      }
    }
  }
}
