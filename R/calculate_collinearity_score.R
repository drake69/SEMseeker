# Function to calculate the collinearity score among covariates
calculate_collinearity_score <- function(df) {

  # Calculate VIF for each variable using a linear model where each variable is regressed on all others
  # Proper VIF calculation
  # vif_values <- tryCatch({
  #   model <- lm( ~ ., data = df)  # Dummy model
  #   car::vif(model)
  # }, error = function(e) {
  #   warning("VIF calculation failed: ", e$message)
  #   return(rep(NA, ncol(df)))
  # })

  # Identify high VIF variables
  # high_vif <- names(vif_values)[vif_values > vif_threshold]

  # vif_values <- apply(df, 2, function(var) {
  #   lm_model <- stats::lm(var ~ ., data = df)
  #   return(max(car::vif(lm_model)))
  # })

  df <- as.data.frame(df)  # Ensure df is a data frame

  # dummy_y <- df[, 1] # Just use the first column as a dummy 'Y' for the lm call.
  # Or, if you prefer:
  dummy_y <- rep(1, nrow(df))

  # Fit a model with ALL your intended predictors and a dummy dependent variable
  # This model will only be used to get the VIF values for the predictors.
  # The formula `dummy_y ~ .` means `dummy_y` regressed on all *other* columns in `predictors_for_vif`
  # It's better to explicitly list the predictors to avoid confusion.
  # `~ .` here refers to all columns in `predictors_for_vif`
  vif_model <- stats::lm(dummy_y ~ ., data = df)

  # Calculate VIF for this model's predictors
  vif_values <- car::vif(vif_model)



  # Your original approach, but with a more direct output of the R-squared
  # This calculates R-squared for each variable when regressed on all others
  r_squared_values <- apply(df, 2, function(var_col) {
    # Create a temporary data frame with the current variable as 'Y'
    # and the rest as 'X' variables
    # find column with the exact values of var_col
    var_col_index <- which(vapply(df, function(x) all(x == var_col), logical(1)))

    temp_df <- data.frame(Y = var_col, df[,-var_col_index])

    # Fit the model: current variable ~ all other variables
    lm_model <- stats::lm(Y ~ ., data = temp_df)

    # Return the R-squared value
    return(summary(lm_model)$r.squared)
  })

  # Now, calculate the VIF based on these R-squared values
  vif_values <- 1 / (1 - r_squared_values)

  removal_values <- colnames(df)[which(vif_values > 10)]  # VIF threshold set to 10


  # Calculate Condition Index
  X <- as.matrix(df)
  condition_index <- kappa(X, exact = TRUE)

  # Perform PCA
  pca <- stats::prcomp(X, scale. = FALSE)
  explained_variance <- summary(pca)$importance[2, ]

  # Normalizing and weighting scores
  vif_score <- mean(vif_values)
  condition_index_score <- condition_index / 30  # Assuming 30 as a threshold
  pca_score <- sum(explained_variance[explained_variance > 0.1])  # Components explaining >10% variance

  # Aggregate score (weights can be adjusted as needed)
  total_score <- (vif_score * 0.33) + (condition_index_score * 0.33) + (pca_score * 0.33)

  if(length(removal_values)!=0) {
    log_event("WARNING: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Collinearity check failed. VIF: ", vif_score, ", Condition Index: ", condition_index_score, ", PCA: ", pca_score)
    log_event("WARNING: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Variables should be removed: ", paste(removal_values, collapse = ", "))
    return(removal_values)
  }
  else
  {
    log_event("DEBUG: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Collinearity check passed. \nVIF: ", vif_score, "\nCondition Index: ", condition_index_score, "\nPCA: ", pca_score)
    log_event("JOURNAL: - Collinearity check passed. \nVIF: ", vif_score, "\nCondition Index: ", condition_index_score, "\nPCA: ", pca_score)
  }
  return(c())
}

