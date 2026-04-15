pow10_model <- function (family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition=samples_sql_condition,key)
{

  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  ssEnv <- get_session_info()
  pow10_params <- unlist(strsplit(as.character(family_test),"_"))
  partition_percentage <- as.numeric(pow10_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  tempDataFrame <- as.data.frame(tempDataFrame)

  dep_var <- sig.formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable
  covariates <- dep_var$covariates

  if(length(pow10_params)==3)
    if(pow10_params[3]=="predictor")
    {
      dependent_variable <- dep_var$independent_variable
      independent_variable <- dep_var$dependent_variable
    }


  # Bias correction to avoid pow10 of 0 or negative numbers
  # bias_dependent_variable <- min(tempDataFrame[, dependent_variable])
  # tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1 - ifelse(bias_dependent_variable < 0, bias_dependent_variable, 0)
  # tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1

  # bias_dependent_variable <- min(tempDataFrame[, dependent_variable])
  # tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1 - ifelse(bias_dependent_variable < 0, bias_dependent_variable, 0)

  # bias_dependent_variable <- min(tempDataFrame[, independent_variable])
  # bias_covariates <- min(tempDataFrame[, covariates])
  # bias <- min(bias_dependent_variable, bias_covariates, 0)
  #
  # if (bias == 0) {
  #   tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1
  #   tempDataFrame[, covariates] <- tempDataFrame[, covariates] + 1
  # }
  # if (bias < 0) {
  #   tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] - bias + 1
  #   tempDataFrame[, covariates] <- tempDataFrame[, covariates] - bias + 1
  # }

  # Split data into training and test sets
  training.samples <- caret::createDataPartition(tempDataFrame[, dependent_variable], p = partition_percentage, list = FALSE)
  train.data <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]


  cov_temp <- covariates
  ### changes
  for (cov in cov_temp) {
    # apply pow10 on independent variable and covariates
    cov_temp[cov_temp==cov] <- paste0("pow(10,", cov, ")")
  }
  ### changes

  # Build dynamic formula for lm (to get starting values)
  if (length(covariates) > 0)
    lm_formula_str <- paste0(dependent_variable, " ~ pow(10,", independent_variable, ") + ", paste(cov_temp, collapse = " + "))
  else
    lm_formula_str <- paste0(dependent_variable, " ~ pow(10,", independent_variable, ")")
  lm_formula <- as.formula(lm_formula_str)
  lm_fit <- stats::lm(lm_formula, data = train.data)

  # Extract starting values for a and b
  start_list <- lm_fit$coefficients
  cov_letters <- letters[1:length(start_list)]
  names(start_list) <- cov_letters

  # Build dynamic formula for nls
  if(length(covariates) == 0)
    rhs <- paste0("a + b * pow(10,", independent_variable, ")")
  else
    rhs <- paste0("a +  pow(10,b * ", independent_variable, ")",paste(paste0(" + pow(10,", cov_letters[2:length(cov_letters)], " * ", covariates, ")"), collapse = ""))
  nls_formula <- as.formula(paste0(dependent_variable, " ~ ", rhs))


  # Fit the NLS model
  log_model_result <- tryCatch({
    stats::nls(
      formula = nls_formula,
      start = start_list,
      data = train.data,
      control = stats::nls.control(warnOnly = TRUE, maxiter = 1000, minFactor = 1e-8)
    )
  }, error = function(e) {
    # Handle error
    return(NA)
  })
  # check if log_model_result is a model
  if (!inherits(log_model_result, "nls"))
    return(res)

  predictions <- stats::predict(log_model_result)
  if(partition_percentage < 1)
  {
    # Make predictions
    predictions_test <- stats::predict(log_model_result, newdata = test.data)
    res <- cbind(res, model_performance(predictions, train.data[,dependent_variable], predictions_test, test.data[,dependent_variable]))
  }
  else
    res <- cbind(res, model_performance(predictions, train.data[,dependent_variable], c(),c()))

  # If model fitting was successful
  tryCatch({
    summary_result <- summary(log_model_result)
    coefficients <- summary_result$coefficients
  }, error = function(e) {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), "Error extracting coefficients from the model: ", e$message)
    return(res)
  })

  if(!exists("summary_result"))
    return(res)

  if (!is.null(summary_result) && ncol(summary_result$coefficients) >= 4) {
    # safe to extract p-values
  } else {
    warning("P-values could not be computed due to singular matrix or failed model summary.")
  }

  rownames(coefficients) <- name_cleaning(c("INTERCEPT", independent_variable, covariates))

  # for a and b extract the p-value
  for (i in 1:nrow(coefficients)) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- name_cleaning(paste0(row_name,"_PVALUE",sep=""))
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name

    pow10_coef_estimate <- data.frame(pow10_a_estimate = coefficients[i,1])
    colnames(pow10_coef_estimate) <- name_cleaning(paste0(row_name,"_ESTIMATE",sep=""))

    res <- cbind(res, p_value, pow10_coef_estimate)
  }

  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  if(plot & !any(is.na(predictions)))
  {

    # predicted <- stats::predict(log_model_result, train.data, type = "response")
    # train.data$predicted <- predicted
    # train.data$fitted <- log_model_result$m$fitted

    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL", name_cleaning(samples_sql_condition)))
    filename  =  file_path_build(chartFolder,
      c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable, covariates, key$COMBINED),
      ssEnv$plot_format)

    ggp <- ggplot2::ggplot(train.data, ggplot2::aes(x =!!ggplot2::sym(independent_variable), y =!!ggplot2::sym(dependent_variable))) +
      ggplot2::geom_point(color = ssEnv$color_palette[1]) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ pow(10,x), se = TRUE, color = ssEnv$color_palette_darker[2]) + # The fitted line will be an exponential curve
      ggplot2::labs(
        title = "", # Puoi lasciare vuoto o aggiungere un titolo
        x = independent_variable,
        y = dependent_variable
      )

    if(partition_percentage < 1)
    {
      ggp <- ggp + ggplot2::geom_point(data = test.data, ggplot2::aes(y = predictions_test, x = eval(parse(text=independent_variable))), color = ssEnv$color_palette_darker[3])+
        ggplot2::geom_point(data = test.data, ggplot2::aes(y = dependent_variable, x = eval(parse(text=independent_variable))), color = ssEnv$ssEnv$color_palette[2])+
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)
    }

    ggplot2::ggsave(
      filename,
      plot = ggp,
      scale = 1,
      width = 9,
      height = 9,
      units = c("in"),
      dpi = as.numeric(ssEnv$plot_resolution_ppi)
    )


  }
  return (res)
}

pow <- function(base, exponent) base^exponent
