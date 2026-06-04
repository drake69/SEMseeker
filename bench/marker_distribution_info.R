# NOTE: currently internal — could be re-exported once @examples are added
marker_distribution_info <- function(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy,
    start_fresh = FALSE, ...)

  # 
  # 

  marker_quantization_metric(result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
  marker_fit_to_theoretical_distribution()

  close_env()

}
