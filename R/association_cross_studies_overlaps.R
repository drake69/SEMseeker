# compare inference associations of differente studies
association_cross_studies_overlaps <- function(inference_detail, studies,alpha = 0.05, adjust_per_area = FALSE,
  adjust_globally = FALSE,pvalue_column="PVALUE_ADJ_ALL_BH",statistic_parameter, adjustment_method = "BH",
  result_folder, ...)
{

  pvalue_column <- name_cleaning(pvalue_column)
  if (nrow(inference_detail) >1)
    inference_detail <- subset(inference_detail, depth_analysis == 3)[1,]

  if (nrow(inference_detail) ==0)
  {
    log_event("ERROR: ",format(Sys.time(), "%a %b %d %X %Y")," No inference details found!")
    return(NULL)
  }

  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  cross_study_env <- ssEnv

  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will perform the cross study association analysys for projects \n ", paste0(studies$STUDY, collapse = ", "))

  color_palette <- ssEnv$color_palette
  localKeys <- unique(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","AREA")])


  aggregated_study_results <- data.frame()
  for (a in 1:nrow(localKeys))
  {
    MARKER <- localKeys[a, "MARKER"]
    AREA <- localKeys[a, "AREA"]
    # for each study in studies
    for (s in 1:nrow(studies))
    {
      # get the inference details for the study
      result_folder_study <- studies[s,"RESULT_FOLDER"]
      ssEnv <- init_env( result_folder =  result_folder_study, start_fresh = FALSE,alpha=alpha, ...)
      temp_res <- association_results_get(inference_detail = inference_detail, marker = MARKER, area= AREA,
        adjust_per_area = adjust_per_area, adjust_globally = adjust_globally, pvalue_column= pvalue_column,
        adjustment_method = adjustment_method, significance = NULL)
      if(nrow(temp_res) != 0)
        temp_res$STUDY <- studies[s,"STUDY"]
      aggregated_study_results <- plyr::rbind.fill(aggregated_study_results, temp_res)
    }
  }

  # change STUDY column to take int account the direction of the statistic
  # aggregated_study_results$STUDY <- as.character(paste0(aggregated_study_results$STUDY,"_", ifelse((aggregated_study_results[,statistic_parameter] > 0),"INCR","DECR") ))
  update_session_info(cross_study_env)
  ssEnv <- get_session_info(result_folder)

  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y")," inference files aggregated!")
  # filename <- inference_file_name( inference_detail, "AGGREGATED",ssEnv$result_folderInference,file_extension = "csv",suffix = "", prefix = "")
  # write.csv2(aggregated_study_results, filename, row.names = F)
  # create a file with max of pvalue column for the cross study and mean of statistic parameter
  if(nrow(aggregated_study_results)>0)
    # write a file for each marker
    for (m in unique(aggregated_study_results$MARKER))
    {
      tt <- aggregated_study_results[aggregated_study_results$MARKER == m, ]
      tt <- tt[tt$DEPTH == 3, ]

      tt$KEY <- paste0(tt$AREA,"_",tt$SUBAREA,"_",tt$MARKER,"_",tt$FIGURE,"_",tt$AREA_OF_TEST)
      SPLIT <- split(tt$KEY, tt$STUDY)
      # get the common keys
      common_keys <- Reduce(intersect, SPLIT)
      # get the common keys
      tt <- tt[tt$KEY %in% common_keys, ]
      # remove KEY column
      tt$KEY <- NULL

      if(statistic_parameter!="")
      {
        tt <- tt[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","DEPTH",statistic_parameter, pvalue_column)]
        # get only "AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","DEPTH" common to STUDY
        tt <- tt %>%
          dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$AREA_OF_TEST, .data$DEPTH) %>%
          dplyr::summarise(
            alpha = max(get(pvalue_column), na.rm = TRUE),
            statistic_parameter = mean(get(statistic_parameter), na.rm = TRUE)
          ) %>%
          dplyr::ungroup()
        colnames(tt)[colnames(tt) == statistic_parameter] <- statistic_parameter
      }
      else
      {
        tt <- tt[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","DEPTH",pvalue_column)]
        # summarise grouping by "AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST" and calculate the max of the pvalues
        tt <- tt %>%
          dplyr::group_by(.data$AREA, .data$SUBAREA, .data$MARKER, .data$FIGURE, .data$AREA_OF_TEST, .data$DEPTH) %>%
          dplyr::summarise(
            alpha = max(get(pvalue_column), na.rm = TRUE)
          ) %>%
          dplyr::ungroup()
      }
      # rename pvalue_column to the original name
      colnames(tt)[colnames(tt) == "alpha"] <- pvalue_column
      dest_folder <- dir_check_and_create(ssEnv$result_folderInference,name_cleaning(inference_detail$areas_sql_condition))
      filename <- inference_file_name(inference_detail, m, dest_folder ,prefix="" )

      if(file.exists(filename))
      {
        old_results <- utils::read.csv2(filename, header = TRUE, stringsAsFactors = FALSE)
        # remove statistic_parameter column
        old_results <- old_results[,!colnames(old_results) %in% c(statistic_parameter)]
        if(!pvalue_column %in% colnames(old_results))
          tt <- merge(old_results, tt, all = TRUE, by = c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","DEPTH"))
        # remove KEY COLUMN
        tt$KEY <- NULL
      }
      utils::write.csv2(tt, filename, row.names = FALSE)
    }

  aggregated_study_results <- data.frame()
  for (signif in c(TRUE,FALSE))
  {
    for (a in 1:nrow(localKeys))
    {
      MARKER <- localKeys[a, "MARKER"]
      AREA <- localKeys[a, "AREA"]
      # for each study in studies
      for (s in 1:nrow(studies))
      {
        # get the inference details for the study
        result_folder_study <- studies[s,"RESULT_FOLDER"]
        ssEnv <- init_env( result_folder =  result_folder_study, start_fresh = FALSE,alpha=alpha, ...)
        temp_res <- association_results_get(inference_detail = inference_detail, marker = MARKER, area= AREA,
          adjust_per_area = adjust_per_area, adjust_globally = adjust_globally, pvalue_column= pvalue_column,
          adjustment_method = adjustment_method, significance = signif)
        if(nrow(temp_res) != 0)
          temp_res$STUDY <- studies[s,"STUDY"]
        aggregated_study_results <- plyr::rbind.fill(aggregated_study_results, temp_res)
      }
    }

    # change STUDY column to take int account the direction of the statistic
    # aggregated_study_results$STUDY <- as.character(paste0(aggregated_study_results$STUDY,"_", ifelse((aggregated_study_results[,statistic_parameter] > 0),"INCR","DECR") ))
    update_session_info(cross_study_env)
    ssEnv <- get_session_info(result_folder)

    if(nrow(aggregated_study_results)==0)
      next

    if(statistic_parameter!="")
      aggregated_study_results <- unique(aggregated_study_results[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STUDY",statistic_parameter, pvalue_column)])
    else
      aggregated_study_results <- unique(aggregated_study_results[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STUDY", pvalue_column)])

    # reshape to a table with study as columns, area as rows and values as cell without aggreagtion
    aggregated_study_results_table <- reshape2::dcast(aggregated_study_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ STUDY, value.var = pvalue_column)
    #
    studies_to_comb <- na.omit(unique(aggregated_study_results$STUDY))
    # calculate the mean of the pvalues
    aggregated_study_results_table[, pvalue_column] <- mean(aggregated_study_results_table[,studies_to_comb], na.rm = TRUE)

    # rename columns
    studies_to_comb_cols <- paste0(studies_to_comb,"_",pvalue_column)
    colnames(aggregated_study_results_table) <- c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST",studies_to_comb_cols, pvalue_column)


    if(statistic_parameter != "")
    {
      aggregated_study_results_table_statistic_parameter <- reshape2::dcast(aggregated_study_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ STUDY, value.var = statistic_parameter)
      aggregated_study_results_table_statistic_parameter[, statistic_parameter] <- mean(aggregated_study_results_table_statistic_parameter[,studies_to_comb], na.rm = TRUE)

      studies_to_comb_cols <- paste0(studies_to_comb,"_",statistic_parameter)
      colnames(aggregated_study_results_table_statistic_parameter) <- c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST",studies_to_comb_cols, statistic_parameter)

      aggregated_study_results_table <- merge(aggregated_study_results_table, aggregated_study_results_table_statistic_parameter, by = c("AREA", "SUBAREA", "MARKER", "FIGURE", "AREA_OF_TEST"))
    }

    aggregated_study_results_table$AREA <- gsub("-", "_", aggregated_study_results_table$AREA)
    aggregated_study_results_table$AREA <- gsub("_", "-", aggregated_study_results_table$AREA)


    markers <- unique(aggregated_study_results_table[, c("MARKER")])
    # for (i in seq_along(markers))
    # {
    #   marker <- markers[i]
    #   filename <- inference_file_name(inference_detail, marker,ssEnv$result_folderInference,file_extension = "csv",suffix = "", prefix = "")
    #   aggregated_study_results_table_marker <- subset(aggregated_study_results_table, MARKER == marker)
    #   aggregated_study_results_table_marker$DEPTH <- 3
    #   write.csv2(aggregated_study_results_table_marker, filename, row.names = F)
    # }
    filename <- inference_file_name(inference_detail, paste0(markers, collapse = "_") ,ssEnv$result_folderInference,file_extension = "csv",suffix = "AGGREGATED", prefix = ifelse(signif,"SIGNIFICANT","NOT_SIGNIFICANT"))
    utils::write.csv2(aggregated_study_results_table, filename, row.names = FALSE)
    log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  aggregated files saved!")
    # for( j in 2:length(studies_to_comb))
    {
      # studies_comb <- combinat::combn(studies_to_comb, j)
      # if (j==length(studies_to_comb))
      #   studies_comb <- data.frame("st"=studies_comb)
      # message(studies_comb)
      # for(k in 1: ncol(studies_comb))
      {


        # results_inference_comb <- subset(aggregated_study_results, STUDY %in% studies_comb[,k])
        results_inference_comb <- aggregated_study_results
        keys <- unique(results_inference_comb[, c("SUBAREA", "AREA", "MARKER", "FIGURE")])
        study_count <- length(na.omit(unique(results_inference_comb$STUDY)))

        # Load library
        # library(VennDiagram)
        for (i in 1:nrow(keys))
        {
          # i <- 1
          area_set <- results_inference_comb[
            results_inference_comb$AREA == keys[i, ]$AREA &
              results_inference_comb$SUBAREA == keys[i, ]$SUBAREA &
              results_inference_comb$MARKER == keys[i, ]$MARKER &
              results_inference_comb$FIGURE == keys[i, ]$FIGURE
            , c("AREA_OF_TEST", "STUDY"), ]
          # table(is.na(results_inference_comb))
          area_set <- na.omit(area_set)
          SPLIT <- split(area_set$AREA_OF_TEST, area_set$STUDY)
          categories <- names(SPLIT)
          if(length(categories)<2)
            next
          # remove _

          categories <- gsub("_", " ", categories)
          folder <- dir_check_and_create(ssEnv$result_folderChart,c("OVERLAPS",areas_sql_condition))
          filename <-
            paste(
              folder,  "/",
              keys[i, ]$AREA,
              "_",
              keys[i, ]$SUBAREA,
              "_",
              keys[i, ]$MARKER,
              "_",
              keys[i, ]$FIGURE,
              "_",
              pvalue_column,
              "_",
              ifelse(signif,"SIGNIFICANT","NOT_SIGNIFICANT"),
              "_",
              alpha,
              "_venn_diagramm.",
              ssEnv$plot_format,
              sep = ""
            )

          # Set up the Venn diagram parameters
          # color_palette <- color_palette[length(SPLIT)]

          # Create a Venn diagram using ggvenn
          venn_plot <- ggvenn::ggvenn(
            data = SPLIT,
            fill_color = color_palette[seq_along(SPLIT)],
            set_name_color = "white"
          )

          # Customize the plot with ggplot2
          venn_plot <- venn_plot +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position = "bottom",
              legend.title = ggplot2::element_blank(),
              legend.text = ggplot2::element_text(size = 10, face = "bold")
            ) +
            ggplot2::scale_fill_identity() +
            ggplot2::scale_fill_manual(
              values = color_palette[seq_along(SPLIT)],
              labels = categories
            ) +
            ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, override.aes = list(color = NA)))

          # Create a Venn diagram using ggvenn
          # venn_plot <- ggvenn::ggvenn(SPLIT,
          #   fill_color = color_palette[seq_along(SPLIT)],
          #   set_name_color = "white")
          #
          # # Customize the plot with ggplot2
          # venn_plot <- venn_plot +
          #   ggplot2::theme_void() +
          #   ggplot2::theme(
          #     legend.position = "bottom",
          #     legend.title = ggplot2::element_blank(),
          #     legend.text = ggplot2::element_text(size = 12, face = "bold")
          #   ) +
          #   ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, override.aes = list(color = NA)))

          # save the plot
          ggplot2::ggsave(filename, venn_plot, width = 2048, height = 2048, dpi = as.numeric(ssEnv$plot_resolution_ppi), unit="px")

          # Chart
          # Set up the Venn diagram parameters
          # color_palette <- color_palette[length(SPLIT)]
          #
          # Plot Venn diagram
          # suppressMessages(
          #   VennDiagram::venn.diagram(
          #     x = SPLIT,
          #     fill = color_palette[seq_along(SPLIT)],
          #     alpha = 0.5,
          #     category.names = categories,
          #     cat.col = rep("black", length(SPLIT)),
          #     cat.cex = 1.2,
          #     cat.fontface = "bold",
          #     cex = 1.5,
          #     output = TRUE,
          #     individuals.in.intersections = TRUE,
          #     disable.logging = TRUE,
          #     filename = filename,
          #     resolution = ssEnv$plot_resolution,
          #     units = "px",
          #     height = 1024,
          #     width = 1024
          #   )
          # )
          log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y"),"  venn diagram completed !")

          overlaps <- Reduce(intersect, SPLIT)
          if(length(overlaps)>0)
          {
            filename <- inference_file_name(inference_detail, paste(keys[i, ]$AREA, keys[i, ]$SUBAREA, keys[i, ]$MARKER, keys[i, ]$FIGURE , sep="_"),ssEnv$result_folderInference,file_extension = "csv",suffix = "OVERLAPS", prefix = ifelse(signif,"SIGNIFICANT","NOT_SIGNIFICANT"))
            overlaps <- data.frame("AREA_OF_TEST" = overlaps)
            utils::write.csv2(overlaps, filename, append = TRUE)
          }
        }
      }
    }

    log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  job completed !")
    # remove all files ending with .log
    files <- list.files(ssEnv$result_folderInference, pattern = ".log")
    for (f in seq_along(files))
    {
      file <- files[f]
      file.remove(paste(ssEnv$result_folderInference, "/", file, sep = ""))
    }
  }
  close_env()
}
