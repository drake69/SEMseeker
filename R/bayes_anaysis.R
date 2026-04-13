#' P_to_be_Case_cond_to_be_Epimutated and sensibility analysis of SEMseeker's mutations and lesions
#'
#' @param result_folder where semseeker's results are stored, the root folder
#' @param maxResources percentage of max system's resource to use
#' @param parallel_strategy which strategy to use for parallel execution see future vignete: possible values, none, multisession,sequential, multicore, cluster
#' @param ... other options to filter elaborations
#'
#' @importFrom doRNG %dorng%
#' @export
bayes_analysis <-
  function(result_folder,
    independent_variable = "Sample_Group",
    maxResources = 90,
    parallel_strategy  = "multicore",
    ...)
  {
    j <- 0
    k <- 0
    z <- 0
    markers <- c("MUTATIONS","LESIONS")
    ssEnv <-
      init_env(
        result_folder =  result_folder,
        maxResources =  maxResources,
        parallel_strategy  =  parallel_strategy,
        start_fresh = FALSE,
        markers = markers,
        ...
      )

    arguments <- list(...)
    areas_selection <- c()
    if (!is.null(arguments[["areas_selection"]]))
    {
      areas_selection <- arguments$areas_selection
    }

    localKeys <- ssEnv$keys_markers_figures
    sample_groups <- c("Reference", "Control", "Case")



    study_summary <-  study_summary_get()
    if (independent_variable=="Sample_Group")
      study_summary <- study_summary[, c("Sample_Group","Sample_ID")]
    else
      study_summary <- study_summary[, c("Sample_Group","Sample_ID",independent_variable)]


    for (a in length(markers))
    {
      # a <- 2
      fileNameResults <- file_path_build(baseFolder =  ssEnv$result_folderEuristic,detailsFilename =  c(markers[a],"bayes_analisys"),extension = "csv")

      localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
      keys <- localKeys_1[localKeys_1$MARKER == markers[a], ]

      # clean keys from already done analysis
      if (file.exists(fileNameResults))
      {
        results <- utils::read.csv2(fileNameResults, header  =  T)
        keys_markers_figures_areas_done <- unlist(apply(unique(results [, c("MARKER", "FIGURE", "AREA", "SUBAREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
        keys_to_be_done <-
          unlist(apply(keys[, c("MARKER", "FIGURE", "AREA", "SUBAREA")], 1, function(x)
            paste(
              x, collapse  =  "_", sep  =  ""
            )))
        keys <-
          keys[!(keys_to_be_done %in% keys_markers_figures_areas_done),]
      }

      nkeys <- nrow(keys)
      variables_to_export_nested <- c("variables_to_export","keys","sample_names","ssEnv","file_path_build")
      if (nrow(keys) > 0)
        for (k in 1:nkeys)
        {
          # k <- 3
          key <- keys [k, ]
          pivot_filename <- pivot_file_name(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)
          if (file.exists(pivot_filename))
          {
            log_event("INFO: ",
              format(Sys.time(), "%a %b %d %X %Y"),
              " Starting to read pivot:",
              pivot_filename,
              ".")
            if(!file.exists(pivot_filename))
              next
            pivot <- readr::read_delim(pivot_file_name,
              col_types = readr::cols(
                .default = readr::col_double(),
                AREA = readr::col_character(),
              ),
              show_col_types=FALSE, progress=FALSE)
            row.names(tempDataFrame) <- tempDataFrame$AREA
            if (length(areas_selection) > 0)
              tempDataFrame <-tempDataFrame[tempDataFrame[, 1] %in% areas_selection,]

            #removes area name (eg. gene...)
            tempDataFrame <- tempDataFrame[, -1]
            if (is.null(dim(tempDataFrame)))
              next
            if (plyr::empty(tempDataFrame) | nrow(tempDataFrame) == 0)
              next
            tempDataFrame <- t(tempDataFrame)
            tempDataFrame <- as.data.frame(tempDataFrame)
            log_event(
              "INFO: ",
              format(Sys.time(), "%a %b %d %X %Y"),
              " Read pivot:",
              pivot_filename,
              " with ",
              ncol(tempDataFrame),
              " rows."
            )
            tempDataFrame$Sample_ID <- rownames(tempDataFrame)
            tempDataFrame <- merge(study_summary,tempDataFrame, by  =  "Sample_ID", all.x=TRUE)
            tempDataFrame <- tempDataFrame[, c(-1,-3)]
            tempDataFrame <-subset(tempDataFrame, "Sample_Group"  !=   "Reference")
            tempDataFrame <-subset(tempDataFrame, "Sample_Group"  !=   0)
            tempDataFrame <- as.data.frame(tempDataFrame)
            tempDataFrame$Sample_Group <- tempDataFrame$Sample_Group =="Case"

            tempDataFrame[is.na(tempDataFrame)] <- 0
            cols <- (gsub(" ", "_", colnames(tempDataFrame)))
            cols <- (gsub("-", "_", cols))
            cols <- (gsub(":", "_", cols))
            cols <- (gsub("/", "_", cols))
            cols <- (gsub("'", "_", cols))
            tempDataFrame <- as.data.frame(tempDataFrame)
            if (length(colnames(tempDataFrame)) !=  length(cols))
              stop("ERROR: I'm stopping here data to analyze are not correct, file a bug!")
            colnames(tempDataFrame) <- cols

            phenotype <- tempDataFrame[,independent_variable]
            n_case <- sum(phenotype)
            n_control <- nrow(tempDataFrame) - n_case

            if(ssEnv$showprogress)
              progress_bar <- progressr::progressor(along = 2:ncol(tempDataFrame))
            else
              progress_bar <- ""

            var_to_export <- c("tempDataFrame","ssEnv","progress_bar","phenotype","keys","k","n_case","n_control")
            # for (c in 2:ncol(tempDataFrame))
            results_temp <- foreach::foreach(c  =  2:ncol(tempDataFrame), .combine  =  rbind, .export  =  var_to_export) %dorng%
              {
                area = names(tempDataFrame)[c]
                if(ssEnv$showprogress)
                {
                  progress_bar(sprintf("genomic area: %s", stringr::str_pad( area, 20, side=c('left'), pad=' ')))
                }
                epimutated <- as.numeric(tempDataFrame[,c])!=0
                # apply bayes theorem
                # P(A|B) = P(B|A) * P(A) / P(B)
                # P(A) = probability to be a case
                # P(B) = probability to have a mutation in a specific area
                # P(B|A) = probability to have a mutation in a specific area given that the sample is a case
                # P(A|B) = probability to be a case given that the sample has a mutation in a specific area
                P_to_be_Epimutated <- sum(epimutated) / length(epimutated)

                P_to_be_Case <- n_case / (n_case + n_control)
                P_to_be_Epimutated_cond_to_be_Case <- sum(epimutated & phenotype) / sum(phenotype)
                P_to_be_Case_cond_to_be_Epimutated <- (P_to_be_Epimutated_cond_to_be_Case * P_to_be_Case) / P_to_be_Epimutated

                P_to_be_Control <- n_control / (n_case + n_control)
                P_to_be_Epimutated_cond_to_be_Control <- sum(epimutated & (!phenotype)) / sum(!phenotype)
                P_to_be_Control_cond_to_be_Epimutated <- (P_to_be_Epimutated_cond_to_be_Control * P_to_be_Control) / P_to_be_Epimutated

                if(P_to_be_Case_cond_to_be_Epimutated >= 0.9 & P_to_be_Control_cond_to_be_Epimutated < 0.1)
                  data.frame(as.character(keys[k, "MARKER"]),as.character(keys[k, "FIGURE"]), as.character(keys[k, "AREA"]), as.character(keys[k, "SUBAREA"]),area, P_to_be_Case_cond_to_be_Epimutated,P_to_be_Control_cond_to_be_Epimutated)
              }
            if (is.null(dim(results_temp)))
              next
            results_temp <- as.data.frame(results_temp)
            colnames(results_temp) <- c("MARKER", "FIGURE", "AREA", "SUBAREA", "AREA_OF_TEST", "P_to_be_Case_cond_to_be_Epimutated","P_to_be_Control_cond_to_be_Epimutated")
            if (exists("results"))
              results <- plyr::rbind.fill(results, results_temp)
            else
            {
              results <- results_temp
            }
            utils::write.csv2(x = results, file = fileNameResults, row.names = FALSE)
          }
        }
      # sort descending by P_to_be_Case_cond_to_be_Epimutated and P_to_be_Case_cond_to_be_Epimutated
      # results <- results[order(-results$P_to_be_Case_cond_to_be_Epimutated, -results$P_to_be_Case_cond_to_be_Epimutated),]
      # drop column if all is na
      if (!exists("results"))
        next
      results <- subset(results, results$AREA != "CHR")
      results <- results[, colSums(is.na(results)) != nrow(results)]
      results$P_to_be_Case_cond_to_be_Epimutated <- as.numeric(results$P_to_be_Case_cond_to_be_Epimutated)
      results$P_to_be_Control_cond_to_be_Epimutated <- as.numeric(results$P_to_be_Control_cond_to_be_Epimutated)
      utils::write.csv2(x = results, file = fileNameResults, row.names = FALSE)

      max_P_to_be_Case_cond_to_be_Epimutated <- max(results$P_to_be_Case_cond_to_be_Epimutated)
      max_P_to_be_Case_cond_to_be_Epimutated <- max(results$P_to_be_Case_cond_to_be_Epimutated)
      results <- subset(results,results$P_to_be_Case_cond_to_be_Epimutated!=0 & results$P_to_be_Case_cond_to_be_Epimutated!=0)
      results <- subset(results,results$P_to_be_Case_cond_to_be_Epimutated== max_P_to_be_Case_cond_to_be_Epimutated)
      fileNameResults <- file_path_build(baseFolder =  ssEnv$result_folderEuristic,detailsFilename =  c(markers[a],"filtered_bayes_analisys"),extension = "csv")
      write.csv2(x = results, file = fileNameResults, row.names = FALSE)

      rm(results)
    }
    close_env()
  }
