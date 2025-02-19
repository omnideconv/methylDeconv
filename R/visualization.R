#' Function to plot the results from `deconvolute()`  as boxplots for each celltype
#'
#' @param result result from `deconvolute()` 
#'
#' @export
#'
#' 
results_boxplot <- function(result){
  
  result |> 
    tibble::rownames_to_column(var = 'sample') |>
    tidyr::pivot_longer(cols = -sample, names_to = 'celltype') |>
    ggplot2::ggplot(mapping = ggplot2::aes(x=value, y= celltype, fill=celltype))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = 'top')
}


#' Function to plot the results from `deconvolute()`  as barplots for each sample
#'
#' @param result result from `deconvolute()` 
#'
#' @export
#'
#' 
results_barplot <- function(result){
  
  result |> 
    tibble::rownames_to_column(var = 'sample') |>
    tidyr::pivot_longer(cols = -sample, names_to = 'celltype') |>
    ggplot2::ggplot(mapping = ggplot2::aes(x=value, y=sample, fill=celltype))+
      ggplot2::geom_col()+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = 'top')
}

#' Function to plot the results from `deconvolute_combined()` as boxplots for each celltype and celltype
#'
#' @param result result from `deconvolute_combined()` 
#'
#' @export
#'
#' 
results_aggregated_boxplot <- function(result){
  
  result |> 
    ggplot2::ggplot(mapping = ggplot2::aes(x=value, y= celltype, fill=method))+
    ggplot2::geom_boxplot()+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = 'top')
}



