phenotype_analysis_name <- function(inference_detail, key,prefix = "", suffix="", pvalue_column, alpha, significance)
{
  ssEnv <- get_session_info()

  covariates <- inference_detail$covariates
  covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  TRUE)))
  covariates_dummy <- inference_detail$covariates_dummy
  covariates_dummy <- if(length(covariates_dummy) !=  0 && !is.null(covariates_dummy)) unlist(t(strsplit( gsub(" ","",covariates_dummy),split  =  "+", fixed  =  TRUE)))
  covariates_pca <- ifelse(is.null(inference_detail$covariates_pca),FALSE,inference_detail$covariates_pca)

  family_test <- inference_detail$family_test
  transformation_y <- inference_detail$transformation_y

  inference_detail$independent_variable <- gsub("_SCALED","",as.character(inference_detail$independent_variable))
  if(inference_detail$transformation_x=="scale")
    inference_detail$independent_variable <- paste0(inference_detail$independent_variable, "_SCALED")
  independent_variable <- split_and_clean(inference_detail$independent_variable)

  significance_label <- ifelse(significance, "significant", "non_significant")

  if(is.null(covariates) || length(covariates)  ==  0)
    file_suffix <- ""
  else
    file_suffix <- paste(covariates, collapse = "_")

  file_suffix = paste(file_suffix, suffix, sep = "_")

  if(length(covariates_dummy)!=0)
    file_suffix <- paste0(file_suffix,paste(covariates_dummy, collapse = "_"), "_dummy", sep = "_")

  if(length(covariates_pca)>0)
    file_suffix <- paste0(file_suffix, "pca", sep = "_")

  analysis_name <- paste(as.character(key$MARKER),as.character(key$FIGURE),as.character(key$AREA),as.character(key$SUBAREA), pvalue_column, significance_label, alpha, prefix , as.character(transformation_y), as.character(family_test), file_suffix, sep="_")

  analysis_name <- name_cleaning(analysis_name)

  return(analysis_name)
}
