# (family_test, tempDataFrame, sig.formula, transformation_y, plot)
assoc_mediation_quantreg_model <- function(family_test,tempDataFrame, sig.formula, transformation_y, plot) {

  ssEnv <- core_get_session_info()

  quantreg_params <- unlist(strsplit(as.character(family_test),"_"))
  tau <- as.numeric(quantreg_params[2])
  res <- data.frame("tau" = tau)
  permutations_test <- as.numeric(quantreg_params[3])
  permutations <- as.numeric(quantreg_params[4])

  # decompose the formula, to move the mediator
  # the expected formula is "mediator ~ outcome + treatment + covariates"
  vars <- assoc_sig_formula_vars(sig.formula)
  mediator <- vars[[1]] # burden_Value
  outcome <- vars[[2]] # dependent variable
  treatment <- vars[[3]][1]

  core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Starting mediation analysis for: ", mediator, " and ", outcome)

  res$mediator <- mediator
  res$outcome <- outcome
  res$treatment <- treatment

  if(assoc_calculate_collinearity_score(df_mediator) > 0.5)
  {
    core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Collinearity detected in the mediator model")
    res$collinearity_mediator <- TRUE
    return(res)
  } else
    res$collinearity_mediator <- FALSE

  if(length(vars[[3]])>1)
  {
    covariates <- vars[[3]][-c(1)]
    # transform covariates in a vector of string
    covariates <- vapply(covariates, as.character, character(1))
  } else
    covariates <- c()

  df_mediator <- data.frame(tempDataFrame)
  df_mediator$mediator <- df_mediator[,mediator]
  df_mediator$treatment <- df_mediator[,treatment]
  if(length(covariates) != 0)
    df_mediator <- df_mediator[,c("treatment", "mediator",covariates)]
  else
    df_mediator <- df_mediator[,c("treatment", "mediator")]

  df_outcome <- data.frame(tempDataFrame)
  df_outcome$outcome <- df_outcome[,outcome]
  df_outcome$mediator <- df_outcome[,mediator]
  df_outcome$treatment <- df_outcome[,treatment]
  if(length(covariates) != 0)
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment",covariates)]
  else
    df_outcome <- df_outcome[,c("outcome", "mediator", "treatment")]

  perms <- sort(unique(c(permutations_test, permutations)))

  for (i in seq_along(perms)) {

    env <- new.env()
    env$mediate_env <- mediation::mediate
    env$perm_env <- perms[i]
    env$tau_env <- tau
    env$df_mediator <- df_mediator
    env$df_outcome <- df_outcome

    mediation_result <- tryCatch({
      R.utils::withTimeout({
        suppressMessages(
          eval(expr = quote(
            mediate_env(
              model.m =  quantreg::rq(formula = as.formula(" mediator ~ . ") , data = df_mediator, tau = 0.5),
              model.y =  quantreg::rq(formula =  as.formula("outcome ~ . ") , data = df_outcome, tau = 0.5),
              treat = "treatment",
              mediator = "mediator",
              boot = TRUE,
              sims = perm_env)),
            envir = env)
        )
      }, timeout = 3, onTimeout = "error")  # Set timeout to 60 seconds
    }, TimeoutException = function(ex) {
      message("DEBUG: Timeout: Mediation analysis took too long.")
      NULL
    }, error = function(e) {
      message("DEBUG: Error in mediation analysis: ", e)
      NULL
    })

    # remove environment
    rm(env)

    if(is.null(mediation_result))
      return(data.frame())

    # suppressMessages(
    #   mediation_result <- eval(expr = quote(
    #     mediate_env(
    #       model.m =  quantreg::rq(formula = as.formula(" mediator ~ . ") , data = df_mediator, tau = 0.5),
    #       model.y =  quantreg::rq(formula =  as.formula("outcome ~ . ") , data = df_outcome, tau = 0.5),
    #       treat = "treatment",
    #       mediator = "mediator",
    #       boot = TRUE,
    #       sims = perm_env)),
    #     envir = env)
    # )


    res$pvalue <- max(mediation_result$d0.p, mediation_result$z0.p, mediation_result$tau.p, mediation_result$n0.p)
    res$permutations <- perms[i]
    # Extract key components from the summary
    res$acme <- mediation_result$d0
    res$acme_ci_start <- mediation_result$d0.ci[1]
    res$acme_ci_end <- mediation_result$d0.ci[2]
    res$acme_pvalue <- mediation_result$d0.p

    res$ade <- mediation_result$z0
    res$ade_ci_start <- mediation_result$z0.ci[1]
    res$ade_ci_end <- mediation_result$z0.ci[2]
    res$ade_pvalue <- mediation_result$z0.p

    res$total_effect <- mediation_result$tau.coef
    res$total_effect_ci_start <- mediation_result$tau.ci[1]
    res$total_effect_ci_end <- mediation_result$tau.ci[2]
    res$total_effect_pvalue <- mediation_result$tau.p

    res$prop_mediated <- mediation_result$n0
    res$prop_mediated_ci_start <- mediation_result$n0.ci[1]
    res$prop_mediated_ci_end <- mediation_result$n0.ci[2]
    res$prop_mediated_pvalue <- mediation_result$n0.p
    if(res$pvalue < as.numeric(ssEnv$alpha))
      next
  }

  core_log_event("DEBUG: ",  format(Sys.time(), "%a %b %d %X %Y"), " - Mediation analysis completed for: ", mediator, " and ", outcome, " with p-value: ", res$pvalue)

  # Return the summary of the mediation analysis
  return(res)
}

