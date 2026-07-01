assoc_sig_formula_vars <- function(sig.formula){

  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]
  # check if contains +
  covariates <- c()
  if(grepl("\\+",independent_variable))
  {
    vars <- strsplit(independent_variable,"\\+")
    vars <- unlist(vars)
    # take the first
    independent_variable <- vars[1]
    # the second to the end
    covariates <- vars[-1]
  }
  return(list(dependent_variable=dependent_variable,independent_variable=independent_variable, covariates=covariates))

}
