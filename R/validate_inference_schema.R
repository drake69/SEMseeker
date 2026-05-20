#' Fill missing optional columns in inference_details with NA
#'
#' Extracted from association_analysis() (was inline at lines 105-113).
#' Ensures downstream code can always reference these columns without
#' having to check for their existence.
#'
#' @param inference_details data.frame from the user.
#' @return data.frame with all expected columns present (NA where missing).
#' @keywords internal
validate_inference_schema <- function(inference_details) {
  inference_details <- as.data.frame(inference_details)
  expected_values <- c("independent_variable", "family_test", "covariates",
    "covariates_dummy", "transformation_y", "depth_analysis", "filter_p_value",
    "samples_sql_condition", "collinearity_check", "covariates_pca")
  missing_cols <- expected_values[!expected_values %in% colnames(inference_details)]
  for (ev in missing_cols) {
    inference_details[, ev] <- NA
  }
  inference_details
}
