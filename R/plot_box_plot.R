plot_box_plot <- function (dataFrameToPlot, independent_variable,dependent_variable, transformation_y, family_test, samples_sql_condition="",key)
{
  if (!assoc_is_family_dicotomic(family_test))
    return()

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- core_get_session_info()
  chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("COMPARISON",core_name_cleaning(as.character(samples_sql_condition))))
  filename  <-  io_file_path_build(chartFolder,toupper(c("BOX_PLOT",family_test,as.character(transformation_y), independent_variable,"Vs", dependent_variable,area, subarea, marker, figure)),ssEnv$plot_format)
  if(!file.exists(filename))
  {
    # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution)
    # formula <- as.formula(paste(dependent_variable,"~",independent_variable, sep=""))
    # dataFrameToPlot[, independent_variable] <- as.factor(dataFrameToPlot[, independent_variable])
    # # dataFrameToPlot$Sample_Group  <- stats::relevel(as.factor(dataFrameToPlot$Sample_Group), "Control")
    # # Number of boxplots
    num_boxplots <- length(unique(dataFrameToPlot[, independent_variable]))
    # # Define the labels for the independent variables
    # # labels <- levels(independent_variable)
    # # replace underscore with space
    # # labels <- gsub("_", " ", labels)
    # # Create the boxplot with specified colors and labels
    # if (num_boxplots > length(ssEnv$color_palette))
    #   graphics::boxplot(formula, data= dataFrameToPlot,   cex = 2)
    # else
    #   graphics::boxplot(formula, data= dataFrameToPlot,   cex = 2, col = ssEnv$color_palette[1:num_boxplots])
    # grDevices::dev.off()

    p <- ggplot2::ggplot(dataFrameToPlot, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = !!ggplot2::sym(independent_variable)), outlier.size = 1.5) +
      ggplot2::coord_cartesian(ylim = c(min(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*0.95,
        max(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*1.05)) +  # Zoom in to remove extremes
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    if (num_boxplots <= length(ssEnv$color_palette)) {
      p <- p + ggplot2::scale_fill_manual(values = ssEnv$color_palette[seq_len(num_boxplots)])
    }

    ## add pvalue
    if (family_test == "t.test") {
      p <- p + ggpubr::stat_compare_means (method = "t.test", label = "p.format", label.y = 0.5, size = 3)
    } else if (family_test == "wilcox.test") {
      p <- p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.5, size = 3)
    } else if (family_test == "kruskal.test") {

      # p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, size = 3)
      # Create the plot
      groups <- unique(dataFrameToPlot[,independent_variable])
      if(length(groups) < 2)
        return()
      pairwise_comparisons <- combn(groups, 2, simplify = FALSE)

      p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, size = 3) +
        ggpubr::stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", size = 3)

    } else if (family_test == "anova") {
      p <- p + ggpubr::stat_compare_means(method = "anova", label = "p.format", label.y = 0.5, size = 3)
    }
    # ggplot2::ggsave(filename = filename, plot = p, width = 2480/ssEnv$plot_resolution, height = 2480/ssEnv$plot_resolution, dpi = as.numeric(ssEnv$plot_resolution_ppi))
    ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, dpi = as.numeric(ssEnv$plot_resolution_ppi), unit="in")


    #### VIOLIN PLOT #####
    filename  <-  io_file_path_build(chartFolder,
      toupper(c("VIOLIN_PLOT",family_test,as.character(transformation_y), independent_variable,"Vs", dependent_variable,area, subarea, marker, figure)),ssEnv$plot_format)

    # p <- ggplot2::ggplot(dataFrameToPlot, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
    #   ggplot2::geom_violin(ggplot2::aes(fill = !!ggplot2::sym(independent_variable))) +
    #   ggplot2::stat_summary(fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1), geom = "crossbar", width = 0.2, color = "black") +
    #   ggplot2::stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
    #   ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +  # Remove blank space
    #   ggplot2::theme_bw(base_size = 15) +
    #   ggplot2::theme(
    #     legend.position = "bottom",
    #     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    #   )

    p <- ggplot2::ggplot(dataFrameToPlot, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_violin(ggplot2::aes(fill = !!ggplot2::sym(independent_variable))) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1), geom = "crossbar", width = 0.2, color = "black") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
      ggplot2::coord_cartesian(ylim = c(min(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*0.95,
        max(dataFrameToPlot[[dependent_variable]], na.rm=TRUE)*1.05)) +  # Zoom in to remove extremes
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )

    if (num_boxplots <= length(ssEnv$color_palette)) {
      p <- p + ggplot2::scale_fill_manual(values = ssEnv$color_palette[seq_len(num_boxplots)])
    }

    ## add pvalue
    if (family_test == "t.test") {
      p <- p + ggpubr::stat_compare_means (method = "t.test", label = "p.format", label.y = 0.5, size = 3)
    } else if (family_test == "wilcox.test") {
      p <- p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.5, size = 3)
    } else if (family_test == "kruskal.test") {

      # p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, size = 3)
      # Create the plot
      groups <- unique(dataFrameToPlot[,independent_variable])
      if(length(groups) < 2)
        return()
      pairwise_comparisons <- combn(groups, 2, simplify = FALSE)

      p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, size = 3) +
        ggpubr::stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", size = 3)

    } else if (family_test == "anova") {
      p <- p + ggpubr::stat_compare_means(method = "anova", label = "p.format", label.y = 0.5, size = 3)
    }
    ggplot2::ggsave(filename = filename, plot = p, width = 7, height = 7, dpi = as.numeric(ssEnv$plot_resolution_ppi), unit="in")

  }
}
