#' Manhattan plot of association results per genomic area
#'
#' @param marker character. SEM metric to plot (e.g. \code{"MUTATIONS"},
#'   \code{"DELTARP"}, \code{"DELTAR"}, \code{"DELTAQ"}).
#' @param figure character. Mutation direction: \code{"HYPO"} or \code{"HYPER"}.
#' @param area character. Genomic area level (e.g. \code{"GENE"},
#'   \code{"ISLAND"}, \code{"DMR"}).
#' @param subarea character. Sub-area within the area (e.g. \code{"TSS1550"},
#'   \code{"WHOLE"}).
#' @param family character. Statistical family used in the association analysis
#'   (e.g. \code{"wilcoxon"}, \code{"gaussian"}).
#' @param adjust_method character. Column name of the adjusted p-value to use
#'   for colouring (e.g. \code{"BH"}).
#' @param phenotype character. Sample sheet column used to colour points
#'   (e.g. \code{"Sample_Group"}).
#' @param only_significant_areas logical. If \code{TRUE}, show only regions
#'   with adjusted p-value < 0.05 (default \code{FALSE}).
#' @return Invisibly \code{NULL}. A Manhattan plot PNG is saved under
#'   \code{Charts/MARKER_PER_AREA/} in the active result folder.
#'
#' @examples
#' result_dir <- tempdir()
#' \dontrun{
#' manhattan_plot_per_area(
#'   marker        = "DELTARP",
#'   figure        = "HYPO",
#'   area          = "GENE",
#'   subarea       = "WHOLE",
#'   family        = "wilcoxon",
#'   adjust_method = "BH",
#'   phenotype     = "Sample_Group"
#' )
#' }
#' @export
#'
manhattan_plot_per_area <- function(marker,figure,area,subarea,family, adjust_method,phenotype, only_significant_areas=FALSE){

  if(only_significant_areas)
  {
    # get inference file

    # select for family, marker and figure, adjust_method, area and subarea
  }

  pivot_data <- pivot_to_long_format(marker, figure, area,subarea)

  group_variable <- pivot_data[,area]

  ggplot2::ggplot(pivot_data, ggplot2::aes(x = factor(group_variable), y = ifelse(pivot_data$VALUE < 0, 0, pivot_data$VALUE), fill = phenotype)) +
    ggplot2::geom_point(shape = 21, position = ggplot2::position_jitter(), alpha = 0.7, color = as.numeric(pivot_data[,area])) +
    ggplot2::labs(x = area, y = "VALUE", title = paste("Outlier per ",area)) +
    ggplot2::scale_fill_gradient(low = "white", high = ssEnv$color_palette[1]) +
    ggplot2::labs(title = paste(marker, " per ", area), x = area, y = marker) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

}
