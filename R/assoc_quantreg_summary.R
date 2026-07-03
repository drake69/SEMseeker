#' Quantile regression result value, confidence interval and p.value
#'
#' @param permutation_vector vector of permutation statistc signal regression
#' @param statistic_parameter signal regression
#' @param conf.level confidence intervals alpha level
#'
#' @return ci and p.value and pvalue_limit
#' @importFrom doRNG %dorng%
assoc_exact_pvalue <-function(permutation_vector, statistic_parameter, conf.level){

  pvalue_limit <- 1 - conf.level
  pvalue_limit_inf <- (pvalue_limit/2)
  pvalue_limit_sup <- 1 - (pvalue_limit/2)

  ci <- stats::quantile(permutation_vector, probs = c(pvalue_limit_inf, 0.5, pvalue_limit_sup))
  ci.lower <- ci[1]
  ci.median <- ci[2]
  ci.upper <- ci[3]
  if (statistic_parameter>ci.median)
    p.value <- mean(permutation_vector>statistic_parameter)
  else
    p.value <- mean(permutation_vector<statistic_parameter)

  return(c(ci.lower,ci.upper, p.value, pvalue_limit))
}
