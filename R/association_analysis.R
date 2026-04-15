#' Association analysis of SEMseeker results
#'
#' Run statistical association models between SEM metrics and a phenotype
#' variable. Supports group tests (Wilcoxon, t-test), GLM families (gaussian,
#' poisson, binomial), quantile regression, correlations (Pearson, Kendall,
#' Spearman), and multi-covariate formulas (e.g.
#' \code{MUTATIONS_* ~ covariate1 + covariate2}).
#'
#' @param inference_details data.frame. Each row defines one analysis run.
#'   Required columns:
#'   \describe{
#'     \item{independent_variable}{Sample sheet column used as grouping /
#'       covariate variable.}
#'     \item{family_test}{Statistical model: \code{"wilcoxon"},
#'       \code{"stats::t.test"}, \code{"gaussian"}, \code{"poisson"},
#'       \code{"binomial"}, \code{"pearson"}, \code{"kendall"},
#'       \code{"spearman"}, or quantile regression as
#'       \code{"quantreg_<tau>_<runs>"} (e.g. \code{"quantreg_0.25_2000"}).}
#'     \item{transformation_y}{Transformation applied to the dependent variable:
#'       \code{"none"}, \code{"scale"}, \code{"log"}, \code{"log2"},
#'       \code{"log10"}, \code{"exp"}, or
#'       \code{"quantile_<n>"} (e.g. \code{"quantile_3"}).}
#'     \item{marker}{SEM metric column prefix (e.g. \code{"DELTARP"},
#'       \code{"MUTATIONS"}).}
#'     \item{depth_analysis}{Integer depth: \code{1} = sample level,
#'       \code{2} = type level (gene, DMR, CpG island),
#'       \code{3} = genomic area (TSS1550, WHOLE, TSS200, …).}
#'   }
#' @param result_folder character. Path to the SEMseeker result folder.
#' @param maxResources numeric. Maximum percentage of CPU cores to use
#'   (default 90).
#' @param parallel_strategy character. Parallelisation backend; possible
#'   values: \code{"none"}, \code{"multisession"}, \code{"sequential"},
#'   \code{"multicore"}, \code{"cluster"} (default \code{"multicore"}).
#' @param start_fresh logical. If \code{TRUE}, delete previous inference
#'   results before running (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{init_env()}.
#'
#' @return Invisibly \code{NULL}. Inference result CSV files are written to
#'   the \code{Inference/} sub-folder of \code{result_folder}, one file per
#'   marker/area/family combination defined in \code{inference_details}.
#' @importFrom doRNG %dorng%
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' association_analysis(
#'   inference_details = data.frame(
#'     independent_variable = "Sample_Group",
#'     family_test          = "wilcoxon",
#'     transformation_y     = "none",
#'     marker               = "DELTARP",
#'     areas                = "GENE"
#'   ),
#'   result_folder     = "~/semseeker_results/",
#'   multiple_test_adj = "BH"
#' )
#' }
#' @export
association_analysis <- function(inference_details,result_folder, maxResources = 90,
  parallel_strategy  = "multicore",start_fresh = FALSE, ...)
{

  j <- 0
  k <- 0
  z <- 0
  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)

  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will perform the association analysys for project \n in ", ssEnv$result_folderData)

  if(start_fresh)
    unlink(ssEnv$result_folderInference, recursive = TRUE)
  dir_check_and_create(ssEnv$result_folderInference,c())

  localKeys <- ssEnv$keys_markers_figures
  sample_groups <- c("Reference","Control","Case")

  deltaX_get()
  annotate_position_pivots()

  inference_details <- unique(inference_details)
  # Assuming your dataset is named 'inference_details'
  cleaned_data <- inference_details %>%
    dplyr::group_by(dplyr::across(-depth_analysis)) %>%  # Group by all columns except 'depth_of_analysis'
    dplyr::slice_max(depth_analysis, n = 1) %>%  # Keep row with max 'depth_of_analysis' per group
    dplyr::ungroup()  # Remove grouping

  # View result
  head(cleaned_data)


  # variables_to_export <- c("n", "working_data", "sig.formula", "tau", "lqm_control", "estimate", "independent_variable", "inference_details", "ssEnv", "%dorng%", "k", "iter", "RNGseed", "checkRNGversion",
  #                          "getRNG", "%||%", ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", ".getRNGattribute", "isNumber", "isReal", "isInteger", ".foreachGlobals", "RNGprovider", ".RNGkind_length", "tail",
  #                          "file_path_build", "dir_check_and_create", "apply_stat_model", "doRNGversion", ".getRNG", "hasRNG", "nextRNG", "RNGkind", "setRNG", "RNGstr")

  # foreach::foreach(z  =  1:nrow(inference_details), .export  =  variables_to_export) %dorng%
  # remove columns with all NA values or empty string
  inference_details <- as.data.frame(inference_details)
  expected_values <- c("independent_variable","family_test","covariates","covariates_dummy","transformation_y",
    "depth_analysis","filter_p_value","samples_sql_condition",
    "collinearity_check","covariates_pca")
  expected_values <- expected_values[!expected_values %in% colnames(inference_details)]
  for (ev in expected_values)
  {
    inference_details[,ev] <- NA
  }

  for(z in 1:nrow(inference_details))
  {
    results <- data.frame()
    start_time <- Sys.time()
    processed_items <- 0
    inference_detail <- inference_details[z,]

    filter_p_value <- if(!is.null(inference_detail$filter_p_value)) inference_detail$filter_p_value else TRUE

    inference_detail_prettified <- t(inference_detail)
    # Generate a plain text table using kable
    inference_detail_prettified <- knitr::kable(inference_detail_prettified, format = "simple",
      align = "l",    # Left align for all columns
      digits = 2,     # Number of digits for numeric columns
      row.names = TRUE) # Suppress row names

    # Join into a single string
    inference_detail_prettified <- paste(inference_detail_prettified, collapse = "\n")
    log_event("JOURNAL: ##############################################################################################################")
    log_event("JOURNAL: ", format(Sys.time(), "%a %b %d %X %Y"), " \nStarting association analysis for inference detail:\n", inference_detail_prettified)

    family_test <- split_and_clean(inference_detail$family_test)

    study_summary <-   study_summary_get(inference_detail$samples_sql_condition)

    results <- data.frame()
    if (validate_family_test(family_test))
    {
      transformation_y <- inference_detail$transformation_y
      if(is.null(transformation_y) || length(transformation_y)  ==  0)
        transformation_y <- NULL

      res_model_covariates <- covariates_model(inference_detail, study_summary)
      study_summary <- res_model_covariates$study_summary
      covariates <- res_model_covariates$covariates
      inference_detail <- res_model_covariates$inference_detail

      independent_variable <- gsub(" ","", inference_detail$independent_variable)

      if(independent_variable %in% covariates)
      {
        log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " The independent variable is also present as covariate!")
        next
      }

      if( is.null(independent_variable) || length(independent_variable)  ==  0)
      {
        log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " One indipendent variable is missed! Skipped.")
      }
      else
      {
        depth_analysis <- inference_detail$depth_analysis
        if( is.null(depth_analysis) || length(depth_analysis)  ==  0)
        {
          depth_analysis <- 1
          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Missed DEPTH analysis inference forced to 1.")
        }

        # transform independent variable as factor
        if(family_test  ==  "binomial" || family_test  ==  "wilcoxon" || family_test  ==  "t.test")
          study_summary[,independent_variable] <- as.factor(study_summary[,independent_variable])

        file_result_prefix <- paste(depth_analysis, as.character(independent_variable),sep = "_")

        #########################################################################################################
        #########################################################################################################

        if (!(independent_variable %in% colnames(study_summary)))
        {
          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " This indipendent variabile:", independent_variable, " is missed! Skipping")
        }
        else
        {
          if(is.null(covariates) || length(covariates)  ==  0)
          {
            sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable)])
          }
          else
          {
            sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable, covariates)])
          }

          if(!is.null(covariates) && !length(covariates)  ==  0)
          {
            sample_names <-   sample_names[, c("Sample_ID", independent_variable, covariates)]
            colnames(sample_names) <- c("Sample_ID", independent_variable, covariates)
          }
          else
          {
            sample_names <-   unique(sample_names[, c("Sample_ID", independent_variable)])
            colnames(sample_names) <- c("Sample_ID", independent_variable)
          }

          # remove samples with missed values
          sample_names <- sample_names[complete.cases(sample_names),]
          if(nrow(sample_names)  ==  0)
          {
            log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " No samples with complete data for the analysis! Skipped.")
            next
          }

          markers <- unique(localKeys$MARKER)
          for (a in seq_along(markers))
          {
            results <- data.frame()
            keys <- localKeys[localKeys$MARKER==markers[a],]
            keys <- unique(keys)
            cols <- keys$COMBINED
            fileNameResults <- inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference,prefix= ifelse(areas_selection==c(),"",paste(areas_selection, "_", sep = "")))
            log_event("JOURNAL:","Result saved into file:", fileNameResults, ".")
            if (sum(cols %in% colnames(study_summary))!=0)
            {
              # temporarily filter to existing columns only
              cols <- cols[cols %in% colnames(study_summary)]
              keys$AREA  =  "SAMPLE_GROUP"
              keys$SUBAREA  =  "SAMPLE"
              g_end <- length(cols)

              if(!is.null(covariates) && !length(covariates)  ==  0)
                study_summary_local <- study_summary[, c(independent_variable, covariates, cols,"Sample_Group") ]
              else
                study_summary_local <- study_summary

              #
              file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size  > 3
              old_results <- data.frame()
              # clean keys from already done association
              if(file_good)
              {
                old_results <- unique(utils::read.csv2(fileNameResults, header  =  TRUE))
                old_results_filtered <- old_results[old_results$DEPTH == 1,]
                keys_markers_figures_areas_done <- unlist(apply(unique(old_results_filtered [, c("MARKER","FIGURE","AREA","SUBAREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys_to_be_done <- unlist(apply(keys[, c("MARKER","FIGURE","AREA","SUBAREA")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys <- keys[!( keys_to_be_done %in% keys_markers_figures_areas_done), ]
              }

              if(nrow(keys)>0)
                for(j in 1:nrow(keys))
                {
                  key  =  keys[j,]
                  key$FIGURE <- as.character(key$FIGURE)
                  key$MARKER <- as.character(key$MARKER)
                  g_start <- 2 + length(covariates)
                  column_selectors <- c(independent_variable, covariates, key$COMBINED)
                  # remove empty items
                  column_selectors <- column_selectors[column_selectors != ""]
                  processed_items <- processed_items + ncol(study_summary_local) - g_start
                  if(any(is.na(study_summary_local[, column_selectors])))
                  {
                    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Missing values in the data frame!")
                    # remove rows with missing values
                    study_summary_local <- study_summary_local[complete.cases(study_summary_local[, column_selectors]),]
                    #
                  }
                  result_temp <- apply_stat_model(tempDataFrame  =  study_summary_local[, column_selectors], g_start  =  g_start , family_test  =  family_test, covariates  =  covariates,
                    key  =  key, transformation_y  =  transformation_y, dototal  =  FALSE, session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,inference_detail$samples_sql_condition, ...)
                  results <- plyr::rbind.fill(results, result_temp)
                }

              if (!is.null(dim(results)))
                if (nrow(results)>0)
                  if("PVALUE_ADJ" %in% colnames(results))
                    results <- results[order(results$PVALUE_ADJ),]

              if (exists("old_results"))
              {
                results <- plyr::rbind.fill(results, old_results)
                rm(old_results)
              }

              results[results == ""] <- NA
              # remove columns where all rows are NA
              results <- results[, colSums(is.na(results)) < nrow(results)]
              results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
              association_analysis_save_results(results,fileNameResults, family_test, filter_p_value )
            }


            # execute for all the areas
            if(depth_analysis>1)
            {
              localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
              keys <- localKeys_1[localKeys_1$MARKER==markers[a],]

              file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size  >10
              dototal <- TRUE

              nkeys <- nrow(keys)
              variables_to_export_nested <- c("variables_to_export", "keys", "sample_names", "independent_variable", "covariates",
                "family_test", "transformation_y", "ssEnv", "depth_analysis","file_path_build", "apply_stat_model")
              if(nrow(keys)>0)
                # result_temp_foreach <- foreach::foreach(k  =  1:nkeys, .combine  =  rbind, .export  =  variables_to_export_nested) %dorng%
                for (k in 1:nkeys)
                {
                  if(exists("tempDataFrame"))
                    rm(list  =  c("tempDataFrame"))
                  key <- keys [k,]
                  if (key$AREA=="POSITION")
                    next
                  pivot_filename <- pivot_file_name_parquet(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)
                  if (file.exists(pivot_filename))
                  {
                    areas_selection_temp <- areas_selection
                    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting to read pivot:", pivot_filename,".")
                    tempDataFrame <- as.data.frame(polars::pl$read_parquet(pivot_filename))

                    if(file.exists(fileNameResults) && file.info(fileNameResults)$size  >10)
                    {
                      old_results <- unique(utils::read.csv2(fileNameResults, header  =  TRUE))
                      area_to_remove <- old_results[old_results$MARKER==key$MARKER & old_results$FIGURE==key$FIGURE
                        & old_results$SUBAREA==key$SUBAREA & old_results$AREA==key$AREA,"AREA_OF_TEST"]
                      tempDataFrame <- tempDataFrame[!(tempDataFrame$AREA %in% area_to_remove),]
                      results <- plyr::rbind.fill(results, old_results)
                      rm(old_results)
                    }
                    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Read pivot:", pivot_filename, " with ", nrow(tempDataFrame), " rows.")
                    tempDataFrame[is.na(tempDataFrame)] <- 0
                    if(length(areas_selection_temp)>0)
                    {
                      #
                      # check if areas_selection_temp is a range
                      if (any(grepl(":",areas_selection_temp)))
                      {
                        areas_selection_temp <- areas_selection_temp[grepl(":",areas_selection_temp)][1]
                        areas_selection_temp <- unlist(strsplit(areas_selection_temp,":"))
                        min_col <- min(as.numeric(areas_selection_temp[1]), nrow(tempDataFrame))
                        max_col <- min(as.numeric(areas_selection_temp[2]), nrow(tempDataFrame))
                        areas_selection_temp <- seq(from  =  min_col, to  = max_col )
                        tempDataFrame <- tempDataFrame[areas_selection_temp, ]
                      }
                      else
                      {
                        # areas_selection_temp <- gsub(":","_",gsub("'","_",gsub("-","_",areas_selection_temp)))
                        # names_area <- gsub(":","_",gsub("'","_",gsub("-","_",tempDataFrame[,1])))
                        tempDataFrame <- tempDataFrame[ tempDataFrame$AREA %in% areas_selection_temp, ]
                      }

                      if(nrow(tempDataFrame) == 0)
                      {
                        log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " No areas selected for the analysis! Skipped.")
                      }
                    }

                    if(is.null(dim(tempDataFrame)))
                      next
                    if(plyr::empty(tempDataFrame) | nrow(tempDataFrame)==0)
                      next

                    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting to execute required test for:", key$MARKER, key$FIGURE, key$AREA, key$SUBAREA,".")

                    batch <- 0
                    chunk_size <- ceiling(6000000/ncol(tempDataFrame))  # Define a chunk size
                    for (i in seq(1, nrow(tempDataFrame), by = chunk_size)) {

                      chunk_indices <- i:min(i + chunk_size - 1, nrow(tempDataFrame))
                      tempDataFrameBatch <- as.data.frame(tempDataFrame)[chunk_indices,]
                      rownames(tempDataFrameBatch) <- tempDataFrameBatch[,1]
                      tempDataFrameBatch <- tempDataFrameBatch[,-1]
                      tempDataFrameBatch <- t(tempDataFrameBatch)
                      tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                      tempDataFrameBatch$Sample_ID <- rownames(tempDataFrameBatch)
                      log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Transposed pivot:", pivot_filename, " with ", ncol(tempDataFrameBatch) -1 , " columns.")
                      batch <- batch + 1
                      if(nrow(tempDataFrameBatch)>1)
                      {
                        tempDataFrameBatch <-  merge( x  =   sample_names, y  =   tempDataFrameBatch,  by.x  =  "Sample_ID",  by.y  =  "Sample_ID" , all.x  =  TRUE)
                        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Merged pivot:", pivot_filename, " with ", ncol(tempDataFrameBatch), " columns.")
                        tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                        tempDataFrameBatch[is.na(tempDataFrameBatch)] <- 0
                        tempDataFrameBatch <- tempDataFrameBatch[,-1]
                        cols <- (colnames(tempDataFrameBatch))
                        tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                        if(length(colnames(tempDataFrameBatch)) !=  length(cols))
                          stop("ERROR: I'm stopping here data to associate are not correct, file a bug!")
                        colnames(tempDataFrameBatch) <- cols
                        # two are the sample_id and the sample_group
                        g_start <- 2 + length(covariates)
                        processed_items <- processed_items + ncol(tempDataFrameBatch) - g_start
                        if(any(is.na(tempDataFrameBatch)))
                        {
                          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Missing values in the data frame!")
                          # next
                        }
                        result_temp_local_batch <- apply_stat_model(tempDataFrame  =  tempDataFrameBatch, g_start  =  g_start, family_test  =  family_test, covariates  =  covariates,
                          key  =  key, transformation_y =  transformation_y, dototal  =  (length(areas_selection_temp)==0),
                          session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,inference_detail$samples_sql_condition, ...)

                        results <- plyr::rbind.fill(results, result_temp_local_batch)
                        results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
                        # save every 100 rows
                      }
                      association_analysis_save_results(results,fileNameResults, family_test, filter_p_value )
                    }

                    association_analysis_log(cbind(inference_detail,keys[k,]), start_time, Sys.time(), processed_items)
                  }
                  else
                    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " File not found:", pivot_filename,".")
                  association_analysis_log(cbind(inference_detail,keys[k,]), start_time, Sys.time(), processed_items)
                  if(nrow(results)!=0)
                    results <- subset(results, MARKER == key$MARKER)
                }

            }

            if(exists("result_temp_foreach"))
              rm(result_temp_foreach)
            if(exists("result_temp_local_batch"))
              rm(result_temp_local_batch)
            if(exists("tempDataFrame"))
              rm(tempDataFrame)

            # remove any column containing in the name samples_sql_condition
            results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
            # results$samples_sql_condition <- inference_detail$samples_sql_condition
            association_analysis_save_results(results,fileNameResults, family_test, filter_p_value )

          }
        }
      }
    }
    if(nrow(results) != 0)
    {
      results$TRANSFORMATION_X <- inference_detail$transformation_x
      association_analysis_save_results(results,fileNameResults, family_test, filter_p_value )
    }
    total_time <- difftime(Sys.time(), start_time, units  =  "mins")
    end_time <- Sys.time()
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Finished processing association analysis for ", processed_items, " items in ", total_time, " minutes." )
    log_event("JOURNAL:", format(Sys.time(), "%a %b %d %X %Y"), " Association Analysis finished in ", total_time, " minutes. \n ####################################################" )
  }
  close_env()
}


