#' Get ggplot Plot Limits
#' 
#' Gets the ranges for the width and height of a ggplot panel.
#' 
#' @param plot A plot created with ggplot2
#'
#' @return A list: xmin, xmax, ymin, ymax giving the coordinates of the limits of a ggplot panel.
#' @details
#' Sometimes when adding text to a ggplot, the text is cutoff if it extends beyond the limits of plot panel. This function provides information enabling the user to extend the panel limits so that the text is not cutoff.
#' 
#' @export
#' @examples
#' library(ggplot2)
#' data(iris)
#' plt <- ggplot(data=iris, aes(x=Species, y=Petal.Length)) + geom_boxplot()
#' get_plot_limits(plt)
get_plot_limits <- function(plot) {
  gb = ggplot2::ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}
