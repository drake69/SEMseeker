#' Title
#'
#' @param family_test test or regression to apply
#' @param transformation_y transformation_y to apply to data
#' @param tempDataFrame data frame to use for test/regression
#' @param independent_variable regressor
#' @param g_start index of the first burden column in tempDataFrame
#' @param g_end index of the last burden column in tempDataFrame
#' @param dototal logical; if TRUE, append a column with the total (row-sum) burden
#' @param covariates vector of covariates to be found in the sample sheet
#' @param depth_analysis 1 only sample, 2 chr, 3 alle genomic areas
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers (used to name TOTAL columns)
#'
#' @return A named list with two elements: \code{tempDataFrame} (the prepared
#'   and optionally transformed data.frame) and \code{independent_variableLevels}
#'   (the factor levels of the independent variable, or \code{NULL} for continuous
#'   outcomes).
#'
data_preparation <- function(family_test,transformation_y,tempDataFrame, independent_variable, g_start, g_end, dototal, covariates, depth_analysis, key)
{

  # browser()
  ssEnv <- get_session_info()

  transformation_y <- as.character(transformation_y)
  independentVariableIsFactor <- FALSE
  independent_variableLevels <- NULL
  if (is.family_dicotomic(family_test))
  {
    test_factor <- as.factor(tempDataFrame[, independent_variable])
    #
    independent_variableLevels <- NA
    independentVariableIsFactor <- FALSE
    if(length(levels(test_factor))>1)
    {
      # sort alphabetically factors
      # levels(test_factor) <- sort(levels(test_factor))
      tempDataFrame[, independent_variable] <- as.factor(tempDataFrame[, independent_variable])
      tempDataFrame[, independent_variable] <- droplevels(tempDataFrame[, independent_variable])
    }
    if(is.factor(tempDataFrame[, independent_variable]))
    {
      independentVariableIsFactor <- TRUE
      independentVariableData <- tempDataFrame[, independent_variable]
      independent_variableLevels <- levels(tempDataFrame[, independent_variable])
      # independent_variable1stLevel <- levels(tempDataFrame[, independent_variable])[1]
      # independent_variable2ndLevel <- levels(tempDataFrame[, independent_variable])[2]
    }
  }
  else
    tempDataFrame <- as.data.frame(sapply(tempDataFrame, as.numeric))

  originalDataFrame <- tempDataFrame
  if (independentVariableIsFactor)
    tempDataFrame[, independent_variable] <- independentVariableData

  df_head <- tempDataFrame[,1:(g_start-1)]

  burden_values <- sapply(tempDataFrame[,g_start:g_end], as.numeric)
  burden_values <- as.data.frame(burden_values)

  df_colnames <- colnames(tempDataFrame)
  if( !is.null(dim(burden_values))  & dototal) {
    sum_area <- apply(burden_values, 1, sum)
    total_label <- paste("TOTAL_",key$MARKER,"_",key$FIGURE, sep="")
    if(depth_analysis==2)
    {
      #select just column of independent variables, remove columns burden value, preserve only total
      df_colnames <- c(df_colnames[!(df_colnames %in% colnames(burden_values))],total_label)
      burden_values <- data.frame(total_label=sum_area)
    }
    else
    {
      burden_values <- data.frame(burden_values,total_label=sum_area)
      df_colnames <- c(df_colnames,total_label)
    }
  }

  if(grepl("log",transformation_y))
  {
    burden_values <- burden_values + min(burden_values[burden_values>0])
  }
  transformation_y <- as.character(transformation_y)
  if(is.null(transformation_y) | length(transformation_y)==0 | is.na(transformation_y))
    transformation_y <- "none"


  burden_values <- as.data.frame(burden_values)
  df_values_orig <- burden_values
  try(
    {
      burden_values = switch(
        as.character(transformation_y),
        "scale" = scale(burden_values),
        "log" = log(burden_values),
        "log2" = log2(burden_values),
        "log10"= log10(burden_values),
        "exp" = exp(burden_values),
        "none" = burden_values,
        burden_values
      )
    }
  )
  if(grepl("quantile", transformation_y))
  {
    qq <- as.numeric(unlist(strsplit(transformation_y,"\\_"))[2])
    burden_values <- as.data.frame(apply(burden_values,2,function(x){
      if(length(unique(x))>=qq)
        as.numeric(dplyr::ntile(x, n=qq))
      else
        rep(0,length(x))
    }))
  }
  burden_values <- as.data.frame(burden_values)

  if(setequal(burden_values,df_values_orig) & transformation_y !="none")
    transformation_y <- paste0("NA_", transformation_y, sep="")


  if(family_test!="binomial" & family_test!="wilcoxon" & family_test!="jsd" & family_test!="t.test" & family_test!="poisson" &
      family_test!="chisq.test" & family_test!="fisher.test" & family_test!="kruskal.test")
  {
    variable_to_transform <- independent_variable
    if(length(covariates)>0)
    {
      variable_to_transform <- c(independent_variable,covariates)
      independent_variableValues <-as.data.frame(apply(tempDataFrame[,variable_to_transform] ,2, as.numeric))
    }
    else
    {
      independent_variableValues <-as.data.frame(as.numeric(tempDataFrame[,variable_to_transform]))
    }
    independent_variableValuesOrig <- independent_variableValues
    suppressWarnings(
      try(
        {
          independent_variableValues = switch(
            as.character(transformation_y),
            "scale" = scale(independent_variableValues),
            "log" = log(independent_variableValues),
            "log2" = log2(independent_variableValues),
            "log10"= log10(independent_variableValues),
            "exp" = exp(independent_variableValues),
            "none" = independent_variableValues,
            independent_variableValues
          )
        }
      )
    )



    # if(grepl("quantile", transformation_y))
    # {
    #   qq <- unlist(strsplit(transformation_y,"_")[2])
    #   df_values_temp <- as.data.frame(apply( burden_values,2,function(x) dplyr::ntile(x, n=qq)))
    #   colnames(df_values_temp) <- colnames(burden_values)
    # }

    if(setequal(burden_values,df_values_orig) & transformation_y !="none")
      transformation_y <- paste0("NA_", transformation_y, sep="")
    else
      tempDataFrame[, variable_to_transform] <- independent_variableValues
  }


  tempDataFrame <- data.frame(df_head, burden_values)
  if(ncol(tempDataFrame)!=length(df_colnames))
    stop("ERROR: I'm stopping here data are not the same size, file a bug!")

  colnames(tempDataFrame) <- df_colnames
  # after the transformation_y some data could be missed
  lost_cols <- colSums(apply(tempDataFrame,2,is.nan))!=0
  lostDataFrame <-  colnames(tempDataFrame)[lost_cols]
  if(sum(lost_cols)!=0)
    utils::write.csv2(lostDataFrame, file.path(ssEnv$session_folder,paste("lost_data_",transformation_y,"_",stringi::stri_rand_strings(1, 12, pattern = "[A-Za-z0-9]"),".log", sep="")))

  #  we want to preserve the NA in the independent variables to be removed by the models
  tempDataFrame[apply(tempDataFrame,2,is.nan)] <- 0
  if(family_test=="binomial")
    tempDataFrame[, independent_variable] <- as.factor(tempDataFrame[, independent_variable])

  # # remove rows with all NA
  # tempDataFrame <- tempDataFrame[,colSums(is.na(tempDataFrame)) != nrow(tempDataFrame)]

  # replace - with _ in colnames
  colnames(tempDataFrame) <- gsub("-", "_", colnames(tempDataFrame))

  result <- list(tempDataFrame, independent_variableLevels)
  names(result) <- c("tempDataFrame", "independent_variableLevels")

  return (result)
}
