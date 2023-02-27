#' Title
#'
#' @param result list of 4
#'
#' @return
#' @export
#'
#' @examples
visualize_result <- function(result){
  if (!is.list(result)) {
    message("this function is for visualizing the this packages result object.")
  }
  if (!is.null(result$epidish)) {
    p1 <- ifelse(is.null(result$epidish$rpc), NULL, visualize_proportions(result$epidish$rpc))
    p2 <- visualize_proportions(result$epidish$cbs)
    p3 <- visualize_proportions(result$epidish$cp)
  }
  if (!is.null(result$tca)) {
    p4 <- visualize_proportions(result$tca)
  }
  if (!is.null(result$flowsortedbloodepic)) {
    p5 <- visualize_proportions(result$flowsortedbloodepic)
  }
  if (!is.null(result$methylcc)) {
    p6 <- visualize_proportions(result$methylcc)
  }
  cowplot::plot_grid(p1, p2, p3, p4, p5, p6)
}

#' Title
#'
#' @param prop Cell Type proportion matrix
#'
#' @return
#' @export
#'
#' @examples
visualize_proportions <- function(prop){
  if (is.data.frame(prop)) {
    return(ggplot2::ggplot(data = reshape::melt(prop),
                           ggplot2::aes(x=variable, y=value))
           + ggplot2::geom_boxplot(ggplot2::aes(fill=variable))
           + ggplot2::theme(legend.position = "none"))
  }
  return(ggplot2::ggplot(data = reshape::melt(prop),
                         ggplot2::aes(x=X2, y=value))
         + ggplot2::geom_boxplot(ggplot2::aes(fill=X2))
         + ggplot2::theme(legend.position = "none"))
}
