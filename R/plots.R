#' Title
#'
#' @param prop Cell Type proportion matrix
#'
#' @return
#' @export
#'
#' @examples
visualize_proportions <- function(prop){
  ggplot2::geom_boxplot(prop)
}
