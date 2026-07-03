assoc_relu_model <- function (family_test, tempDataFrame, sig.formula, transformation_y, plot, samples_sql_condition)
{

  #
  ssEnv <- core_get_session_info()
  relu_params <- unlist(strsplit(as.character(family_test),"_"))
  # core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " exp model fitting executing: ", as.character(sig.formula))
  partition_percentage <- as.numeric(relu_params[2])

  res<- data.frame("PARTITION_PERC" = as.numeric(partition_percentage))

  tempDataFrame <- as.data.frame(tempDataFrame)

  dep_var <- assoc_sig_formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable

  if(length(relu_params)==3)
    if(relu_params[3]=="predictor")
    {
      dependent_variable <- dep_var$independent_variable
      independent_variable <- dep_var$dependent_variable
    }

  # Define the ReLU function
  relu <- function(x) {
    return(pmax(0, x))
  }

  bias_dependent_variable <- min(tempDataFrame[, dependent_variable])
  tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1 - ifelse(bias_dependent_variable < 0, bias_dependent_variable, 0)
  # tempDataFrame[, dependent_variable] <- tempDataFrame[, dependent_variable] + 1

  # bias_independent_variable <- min(tempDataFrame[, independent_variable])
  # tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1 - ifelse(bias_independent_variable < 0, bias_independent_variable, 0)
  tempDataFrame[, independent_variable] <- tempDataFrame[, independent_variable] + 1

  # Split the data into training and test set
  training.samples <- caret::createDataPartition(tempDataFrame[, dependent_variable], p=partition_percentage, list=FALSE)

  train.data  <- tempDataFrame[training.samples, ]
  test.data <- tempDataFrame[-training.samples, ]

  # Assuming you have vectors train.data[,independent_variable] and train.data[,dependent_variable]
  start_a <- exp(stats::coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[1])  # Starting value for a, based on linear fit
  start_b <- coef(stats::lm(log(train.data[,dependent_variable]) ~ train.data[,independent_variable]))[2]       # Starting value for b, based on linear fit

  relu_model_result <- NA
  tryCatch({
    # # Fit the relu model
    formula <- as.formula(paste0(dependent_variable,"relu(",independent_variable,")"))
    relu_model_result <- lm(formula = formula)

  }, error = function(e) {

  })

  predictions <- stats::predict(relu_model_result)
  if(partition_percentage < 1)
  {
    # Make predictions
    predictions_test <- stats::predict(relu_model_result, newdata = test.data)
    res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], predictions_test, test.data[,dependent_variable]))
  }
  else
    res <- cbind(res, assoc_model_performance(predictions, train.data[,dependent_variable], c(),c()))

  # Coefficients and Confidence Intervals
  coefficients <- coef(summary(relu_model_result))
  # conf_int <- confint(relu_model_result)

  # for a and b extract the p-value
  for (i in seq_len(nrow(coefficients))) {
    # i <- 1
    p_value <- coefficients[i,4]
    row_name <- rownames(coefficients)[i]
    pval_name <- paste0("relu_",row_name,"_PVALUE",sep="")
    p_value <- data.frame(p_value)
    colnames(p_value) <- pval_name

    relu_coef_estimate <- data.frame(relu_a_estimate = coefficients[i,1])
    colnames(relu_coef_estimate) <- toupper(paste0("relu_",row_name,"_ESTIMATE",sep=""))

    res <- cbind(res, p_value, relu_coef_estimate)
  }

  # remove rowname from res
  rownames(res) <- NULL
  # res$sig.fromula <- as.character(sig.formula)

  if(plot & !any(is.na(predictions)))
  {
    chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL", core_name_cleaning(samples_sql_condition)))
    filename  <-  io_file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable),ssEnv$plot_format)

    ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
      ggplot2::geom_point(color = ssEnv$color_palette[1]) +
      ggplot2::stat_function(na.rm = TRUE, fun = function(x) (res$relu_A_ESTIMATE * exp(res$relu_B_ESTIMATE * x)), color = ssEnv$color_palette_darker[3]) +
      ggplot2::stat_smooth(method = "loess", color = ssEnv$color_palette[2], se = TRUE) +
      ggplot2::ggtitle("") +
      ggplot2::xlab(independent_variable) +
      ggplot2::ylab(dependent_variable)

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
      width = 1240,
      height = 1240,
      units = c("px"),
      dpi = as.numeric(ssEnv$plot_resolution_ppi)
    )


  }
  return (res)
}

# relu_model_result("relu_0.8", sample_sheet, as.formula("DELTARQ_HYPO ~ Age"))
