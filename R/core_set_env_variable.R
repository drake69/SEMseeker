core_set_env_variable <- function(arguments, var_name, var_default_value, check_values=NULL)
{


  ssEnv <- core_get_session_info()
  if(is.null(ssEnv[[var_name]]))
    ssEnv[[var_name]]  <- var_default_value

  if(!is.null(arguments[[var_name]]))
  {
    if(any(arguments[[var_name]]=="DEFAULT"))
      ssEnv[[var_name]]  <- var_default_value
    else
      ssEnv[[var_name]]  <- arguments[[var_name]]
  }
  core_update_session_info(ssEnv)
  core_log_event("INFO:",format(Sys.time(), "%a %b %d %X %Y"), " " , var_name ," ",  ssEnv[[var_name]])

  if(!is.null(check_values))
  {
    values <- arguments[[var_name]]
    unknown_options <- values[!values %in% check_values]
  }
  else
  {
    unknown_options <- NULL
  }

  # remove from arguments var_name
  arguments[[var_name]] <- unknown_options

  return(arguments)
}
