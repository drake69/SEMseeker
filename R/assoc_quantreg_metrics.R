# extract quantgile regression metrics and plot the results of the quantile regression and the permutation vector histogram
assoc_quantreg_metrics <- function(predicted_values, expected_values, tau, res, family_test,independent_variable,
  transformation_y, dependent_variable, permutation_vector=c(), plot = FALSE ){

  ###########

  ssEnv <- core_get_session_info()

  pinball_loss <- function(expected_values, predicted_values, tau) {
    residuals <- expected_values - predicted_values
    loss <- sum(ifelse(residuals < 0, (1 - tau) * (-residuals), tau * residuals))
    return(loss)
  }

  # Calculate the Pinball Loss
  local_res <- data.frame("pinball_loss" = pinball_loss(expected_values, predicted_values, tau))

  ###########

  # Calculate the proportion of observations below the predicted quantile
  local_res$below_quantile <- mean(expected_values <= predicted_values)

  local_res$qq_distance <-  abs(tau - local_res$below_quantile)

  ###########

  # if(plot & 1!=1)
  # {
  #   chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
  #   filename  =  io_file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable),ssEnv$plot_format)
  #
  #   # Save the plot
  #   # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution)
  #   if(ssEnv$plot_format == "png")
  #     grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution, bg = "transparent")
  #   if(ssEnv$plot_format == "eps")
  #     grDevices::postscript(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution, bg = "transparent")
  #   graphics::plot(expected_values, predicted_values, main = "Observed vs. Predicted",
  #     xlab = "Observed", ylab = "Predicted", pch = 19)
  #   graphics::abline(0, 1, col = ssEnv$color_palette_darker[1] )  # Adding a 45-degree line for reference
  #   grDevices::dev.off()
  #
  #   # plot histograms of betas
  #   if (length(permutation_vector) > 0)
  #   {
  #     chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
  #     filename  =  io_file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation_y), dependent_variable, "HISTOGRAM"),ssEnv$plot_format)
  #
  #     # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution)
  #     if(ssEnv$plot_format == "png")
  #       grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution, bg = "transparent")
  #     if(ssEnv$plot_format == "eps")
  #       grDevices::postscript(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution, bg = "transparent")
  #     graphics::hist(permutation_vector, main = "Histogram of Betas",
  #       xlab = "Regression Coefficent Value", ylab = "Frequency", col = ssEnv$color_palette[1])
  #     grDevices::dev.off()
  #   }
  # }

  if (plot) {
    if (!exists("expected_values") || !exists("predicted_values") || !exists("permutation_vector")) {
      stop("Required variables (expected_values, predicted_values, permutation_vector) are missing.")
    }

    if (length(expected_values) != length(predicted_values)) {
      stop("Length of expected_values and predicted_values must be the same.")
    }

    if (length(expected_values) == 0 || length(predicted_values) == 0) {
      warning("No data available for plotting observed vs. predicted.")
    } else {
      chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart, c("FITTED_MODEL"))
      filename <- io_file_path_build(chartFolder, c(as.character(family_test), independent_variable, "Vs", as.character(transformation_y), dependent_variable), ssEnv$plot_format)

      # Observed vs. Predicted plot
      observed_predicted_plot <- ggplot2::ggplot(data = data.frame(expected_values, predicted_values), ggplot2::aes(x = expected_values, y = predicted_values)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = ssEnv$color_palette_darker[1]) +
        ggplot2::labs(title = "Observed vs. Predicted", x = "Observed", y = "Predicted") +
        ggplot2::theme_minimal()

      ggplot2::ggsave(filename = filename, plot = observed_predicted_plot, device = ssEnv$plot_format, width = 8.27, height = 8.27, dpi = as.numeric(ssEnv$plot_resolution_ppi), bg = "transparent")
    }

    # Check if the permutation_vector has data before plotting histogram
    if (length(permutation_vector) == 0) {
      warning("No data available for plotting histogram of betas.")
    } else {
      chartFolder <- io_dir_check_and_create(ssEnv$result_folderChart, c("FITTED_MODEL"))
      filename <- io_file_path_build(chartFolder, c(as.character(family_test), independent_variable, "Vs", as.character(transformation_y), dependent_variable, "HISTOGRAM"), ssEnv$plot_format)

      betas_histogram <- ggplot2::ggplot(data = data.frame(permutation_vector), ggplot2::aes(x = permutation_vector)) +
        ggplot2::geom_histogram(fill = ssEnv$color_palette[1], color = "black", bins = 30) +
        ggplot2::labs(title = "Histogram of Betas", x = "Regression Coefficient Value", y = "Frequency") +
        ggplot2::theme_minimal()

      ggplot2::ggsave(filename = filename, plot = betas_histogram, device = ssEnv$plot_format, width = 8.27, height = 8.27, dpi = as.numeric(ssEnv$plot_resolution_ppi), bg = "transparent")
    }
  }

  ###########

  # Calculate correlation
  correlation <- stats::cor(expected_values, predicted_values)
  local_res$qq_correlation <- correlation

  ###########
  model_performance_res <- assoc_model_performance(predicted_values, expected_values, c(),c())
  res <- cbind(res, local_res, model_performance_res)

  return(res)
}
