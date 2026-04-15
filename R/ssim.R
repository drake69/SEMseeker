# Compute the Structural Similarity Index (SSIM) between two vectors
ssim <- function(x, y, K1 = 0.01, K2 = 0.03, L = 1) {

  
  are_matrix <- is.matrix(x) & is.matrix(y)
  if(are_matrix)
    # Ensure both vectors have the same length
    if (nrow(x) != nrow(y) | ncol(x)!=ncol(y)) {
      stop("I vettori devono avere la stessa lunghezza.")
    }
  else
    if (length(x) != length(y)) {
      stop("I vettori devono avere la stessa lunghezza.")
    }
  # Compute the means
  mu_x <- mean(x)
  mu_y <- mean(y)


  if(are_matrix)
  {
    n_item <- nrow(x)*ncol(x)
    # Calculate the variance using the formula
    sigma_x2 <- sum((x - mu_x)^2) / ( n_item - 1)
    sigma_y2 <- sum((y - mu_y)^2) / (n_item - 1)
    # Calculate the covariance using the formula
    sigma_xy <- sum((x - mu_x) * (y - mu_y)) / (n_item - 1)
  }
  else
  {
    # Compute the variances
    sigma_x2 <- var(x)
    sigma_y2 <- var(y)
    # Compute the covariance
    sigma_xy <- cov(x, y)
  }



  # Constants to stabilise division with small denominators
  C1 <- (K1 * L)^2
  C2 <- (K2 * L)^2

  # Compute the three SSIM components
  luminance <- (2 * mu_x * mu_y + C1) / (mu_x^2 + mu_y^2 + C1)
  contrast <- (2 * sqrt(sigma_x2) * sqrt(sigma_y2) + C2) / (sigma_x2 + sigma_y2 + C2)
  structure <- (sigma_xy + C2 / 2) / (sqrt(sigma_x2) * sqrt(sigma_y2) + C2 / 2)

  # Compute final SSIM value
  ssim_value <- luminance * contrast * structure
  return(ssim_value)
}
