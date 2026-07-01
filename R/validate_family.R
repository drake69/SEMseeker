validate_family_test <- function(family_test){

  if( is.null(family_test) || length(family_test)  ==  0 || is.na(family_test))
  {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " One test family_test is missed! Skipped.", family_test)
    return(FALSE)
  }

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " family_test: " , as.character(family_test))


  if(family_test=="multinomial" | family_test=="binomial" | family_test=="binomial_bulk" | family_test=="wilcoxon" | family_test=="jsd" | family_test=="t.test" | family_test=="poisson" |
      family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test" | family_test=="pearson" | family_test=="kendall" | family_test=="spearman" |
      family_test=="wilcoxon" | family_test=="gaussian")
    return(TRUE)

  if(grepl("mean-permutation",family_test))
  {
    mean_params <- unlist(strsplit(as.character(family_test),"_"))
    if (length(mean_params) != 4)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "mean-permutation family_test must have been with the follwing syntax mean-permutation_n.permutations.test_n.permutations_conf.level")
    else
      return(TRUE)
  }

  if(grepl("wilcoxon.paired",family_test))
  {
    mean_params <- unlist(strsplit(as.character(family_test),"@"))
    if (length(mean_params) != 2)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "wilcoxon.paired family_test must have been with the following syntax wilcoxon.paired@pairing_variable")
    else
      return(TRUE)
  }

  if(grepl("t.test.paired",family_test))
  {
    mean_params <- unlist(strsplit(as.character(family_test),"@"))
    if (length(mean_params) != 2)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "t.test.paired family_test must have been with the following syntax t.test.paired@pairing_variable")
    else
      return(TRUE)
  }

  if (grepl("quantile-permutation",family_test))
  {
    quantile_params <- unlist(strsplit(as.character(family_test),"_"))
    if (length(quantile_params) == 5)
      return(TRUE)
    core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "quantile-permutation family_test must have been with the follwing syntax quantile-permutation_quantile_n.permutations.test_n.permutations_conf.level")
    return(FALSE)
  }

  if (grepl("quantreg-permutation",family_test))
  {
    quantile_params <- unlist(strsplit(as.character(family_test),"_"))
    if (length(quantile_params) != 5)
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "quantreg-permutation family_test must have been with the follwing syntax quantile-permutation_tau_n.permutations.test_n.permutations_conf.level")
    else
      return(TRUE)
  }

  if (grepl("mediation-ridge", family_test))
    return(TRUE)

  if (grepl("mediation-linear", family_test))
    return(TRUE)

  if (grepl("spearman-permutation",family_test))
    return(TRUE)

  if (grepl("quantreg", family_test))
    return(TRUE)

  if(grepl("polynomial_",family_test))
    return(TRUE)

  # AI-040: limma_<degree>[_<partition>] and voom_<degree>[_<partition>].
  # Same parser shape as polynomial; the actual guard against missing
  # limma installation lives at the dispatch point (apply_stat_model
  # for the batch path, execute_model for the per-area path) following
  # the AI-038 dispatch=guard convention.
  if (grepl("^(limma|voom)_", family_test))
    return(TRUE)

  if (grepl("exp_",family_test))
    return(TRUE)

  if (grepl("log_",family_test))
    return(TRUE)

  if (grepl("log10_",family_test))
    return(TRUE)

  if (grepl("pow10_",family_test))
    return(TRUE)

  core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), "family_test is not recognized: ", family_test)
  return(FALSE)
}
