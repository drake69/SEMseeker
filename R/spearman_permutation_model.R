#' Compute Spearman permutation statistic (CPU)
#'
#' @param sig.formula formula to apply
#' @param df dataframe to use
#' @param shuffle logical; if TRUE, permute the independent variable before computing the correlation
#'
#' @return A numeric scalar: the Spearman correlation coefficient between the
#'   burden and the independent variable (used as the permutation test statistic).
#'
compute_spearman_permutation <- function(sig.formula,df, shuffle = FALSE)
{
  cols <- colnames(df)
  indepVar <- as.character(all.vars(sig.formula)[2])
  burden <- as.character(all.vars(sig.formula)[1])
  tempDataFrame <- df
  idx = which(colnames(tempDataFrame) == indepVar)
  tempDataFrame[,idx] <- as.factor(tempDataFrame[,idx])
  if (shuffle == TRUE)
    tempDataFrame[ ,indepVar] <- sample(tempDataFrame[,indepVar])
  tempDataFrame <- as.data.frame(tempDataFrame)
  colnames(tempDataFrame) <- cols
  # calculate spearman correlation between indepVar and burden
  spearman_coeff <- stats::cor.test(as.numeric(tempDataFrame[,burden]), as.numeric(tempDataFrame[,indepVar]), method="spearman")
  statistic_parameter <- spearman_coeff$estimate
  #
  statistic_parameter <- as.numeric(statistic_parameter[[1]])
  return(statistic_parameter)
}



#' Spearman permutation model
#'
#' @param family_test family test string encoding model type and permutation parameters
#' @param sig.formula formula of the model
#' @param tempDataFrame data
#' @param independent_variable name of regressor
#' @param plot logical; if TRUE, generate diagnostic plots
#' @param samples_sql_condition SQL condition string used to filter samples (used for plot file naming)
#' @param key named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this test
#'
#' @return A numeric p-value from the permutation-based Spearman correlation test.
#'
spearman_permutation <- function(family_test, sig.formula, tempDataFrame, independent_variable,plot, samples_sql_condition=samples_sql_condition, key)
{
  area <- as.character(key$AREA)
  subarea <- as.character(key$SUBAREA)
  marker <- as.character(key$MARKER)
  figure <- as.character(key$FIGURE)

  spearman_params <- unlist(strsplit(as.character(family_test),"_"))
  res <- data.frame()
  # spearman_params template spearman + first_round_of_permutations + second_round_of_permutations + confidence_interval_of_beta
  # apply permutation to obtain signal using spearman
  # Define function to compute delta spearman regression coefficient
  n_permutations_test <- as.numeric(spearman_params[2])
  n_permutations <- as.numeric(spearman_params[3])
  conf.level <- as.numeric(spearman_params[4])
  res$conf.level <- conf.level

  pvalue_limit <- 1 - conf.level
  pvalue_limit_inf <- (pvalue_limit/2)
  pvalue_limit_sup <- 1 - (pvalue_limit/2)

  # Compute signal and p-value for n_permutations replications
  statistic_parameter <-  compute_spearman_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle = FALSE)
  res$statistic_parameter <- statistic_parameter

  permutation_vector <- replicate(n_permutations_test, compute_spearman_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle = TRUE))

  pvalue <- mean(abs(permutation_vector) >= abs(statistic_parameter))
  if (pvalue>1)
    pvalue <- 1
  # Compute average signal and p-value
  if ((pvalue < as.numeric(ssEnv$alpha)) && (n_permutations_test < n_permutations))
    permutation_vector <- replicate(n_permutations, compute_spearman_permutation(sig.formula=sig.formula, df=tempDataFrame, shuffle=TRUE))
  res$r_model <- "spearman_permutation"

  if(length(permutation_vector) == n_permutations_test)
    n_permutations <- n_permutations_test

  res$n_permutations <- n_permutations
  pvalue <- mean(abs(permutation_vector) >= abs(statistic_parameter))
  res$pvalue <- pvalue
  ci <- stats::quantile(permutation_vector, probs = c(pvalue_limit_inf, pvalue_limit_sup))
  ci.lower <- ci[1]
  res$ci.lower <- ci.lower
  ci.upper <- ci[2]
  res$ci.upper <- ci.upper

  return (res)
}
