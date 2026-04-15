#' Filter metrics by transformation type
#'
#' Removes scale-sensitive metrics from the requested set when a non-scale
#' transformation (e.g. \code{"log"}, \code{"sqrt"}) is applied to the
#' dependent variable.  Scale-sensitive metrics (e.g. MAE, RMSE) are
#' meaningless after a non-linear transformation because their units change.
#'
#' @param metrics Character vector of metric names to filter (upper-case).
#' @param transformation_y Character scalar describing the transformation
#'   applied to the dependent variable.  Use \code{"none"} to return all
#'   metrics unchanged, \code{"scale"} to keep all metrics (z-score does not
#'   change units), or any other value (e.g. \code{"log"}) to drop
#'   scale-affected metrics.
#'
#' @return A sorted character vector of metric names that are valid for the
#'   given transformation.
#'
#' @examples
#' SEMseeker:::metrics_filter(c("MAE", "RMSE", "COUNT_SIGN"), "none")
#' SEMseeker:::metrics_filter(c("MAE", "RMSE", "COUNT_SIGN"), "log")
#'
metrics_filter <- function(metrics, transformation_y){

  if(transformation_y=="none")
    return(metrics)

  affected_by_tranformation <- toupper(SEMseeker::metrics_properties[SEMseeker::metrics_properties$Affected_by_Scaling==TRUE,"Metric"])

  if (transformation_y!="scale")
  {
    to_remove_metrics <- metrics %in% affected_by_tranformation
    metrics <- metrics[!to_remove_metrics]
  }
  metrics <- sort(unique(c(metrics)))

  return(metrics)

}
