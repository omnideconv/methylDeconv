#' Function to plot the results from `deconvolute()` or `run_all_methods()` as boxplots.
#'
#' @param result result from `deconvolute()` or `run_all_methods()`
#'
#' @return
#' @export
#'
#' @examples
#' 
#' mset <- minfiData::MsetEx
#' 
#' # one result
#' res_epidish <- deconvolute(mSet, "epidish")
#' visualize_results(res_epidish)
#' 
#' # all results
#' results <- run_all_methods(mSet)
#' visualize_results(res_epidish)
#' 
#' # combine mutliple results into a list
#' res_methylcc <- deconvolute(mset, "methylcc")
#' visualize_results(list(epidish = res_epidish, 
#'                        methycc = res_methylcc))
#' 
visualize_results <- function(result){
  if (is.matrix(result)){
    df <- as.data.frame(result)
    ggplot2::ggplot(data = reshape::melt(df), ggplot2::aes(x=variable,y=value)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill=variable)) +
      ggplot2::labs(x = "Cell Type", y = "proportions") +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90))
  } else if (is.data.frame(result)){
    ggplot2::ggplot(data = reshape::melt(result, id.vars = c("method", "sample")), 
                    ggplot2::aes(x = variable, y = value)) +
     ggplot2::geom_boxplot(ggplot2::aes(fill=variable)) +
     ggplot2::facet_wrap(~method) +
     ggplot2::labs(x = "Cell Type", y = "proportions") +
     ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90))
  } else if (is.list(result)){
    df <- data.frame()
    sapply(names(result), function(m){
      df_tmp <- as.data.frame(result[[m]])
      df_tmp <- rename_cell_types(df_tmp)
      df_tmp$other <- NULL
      df_tmp$method <- m
      df_tmp$sample <- rownames(df_tmp)
      df <<- dplyr::bind_rows(df, df_tmp)
    })

    ggplot2::ggplot(data = reshape::melt(df, id.vars = c("method", "sample")), 
                    ggplot2::aes(x = variable, y = value)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill=variable)) +
      ggplot2::facet_wrap(~method) +
      ggplot2::labs(x = "Cell Type", y = "proportions") +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90))
    
  } else {
    message("this function is for visualizing the this packages result object.")
    return()
  }
}


#' Function to visualize specific cell types of result from `deconvolute()` as box plot.
#'
#' @param result object from `deconvolute()`
#' @param CT vector of cell types to plot, have to be valid column names of `result`
#'
#' @return
#' @export
#'
#' @examples
#' 
#' mset <- minfiData::MsetEx
#'
#' res_epidish <- deconvolute(mSet, "epidish")
#' visualize_result_box(res_epidish, CT = c("B", "CD4T"))
#' 
visualize_result_box <- function(result, CT='B') {
  if (!is.matrix(result)) {
    message("this function is for visualizing one of the this packages result data frames.")
    return()
  }
  df <- as.data.frame(result[,CT])
  
  df <- reshape::melt(df)
  colnames(df) <- c("CT", "value")
  ggplot2::ggplot(data = df, ggplot2::aes(x=CT,y=value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=CT)) +
    ggplot2::labs(y = "proportions") +
    ggplot2::theme(legend.position = "none")
}


#' Function to visualize the results from `deconvolute()` as a stacked bar plot.
#'
#' @param result result object from `deconvolute()`
#' @param rename if the cell types should be filtered and renamed for better comparison (default: FALSE).
#' Might be most useful for MethAtlas. 
#'
#' @return
#' @export
#'
#' @examples
#' mset <- minfiData::MsetEx
#' 
#' res_epidish <- deconvolute(mSet, "epidish")
#' visualize_result_bar(res_epidish)
#' visualize_result_bar(res_methylcc, rename = TRUE)
visualize_result_bar <- function(result, rename = F) {
  if (!is.matrix(result)) {
    message("this function is for visualizing one of the this packages result data frames.")
    return()
  }
  df <- as.data.frame(result)
  if(rename){
    df <- rename_cell_types(df)
  }

  df$sample <- rownames(df)
  ggplot2::ggplot(reshape::melt(df), ggplot2::aes(x = sample, y = value, fill = variable)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::xlab("sample") +
    ggplot2::ylab("Cell Type probability distribution") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  
}


#' Show scatterplot that compares two results.
#'
#' @param res1 result1 from `deconvolute()`
#' @param res2 result2 from `deconvolute()`
#'
#' @return scatter plots, one for each common cell type
#' @export
#'
#' @examples
#' mset <- minfiData::MsetEx
#'
#' res_epidish <- deconvolute(mset, "epidish")
#' res_methylcc <- deconvolute(mset, "methylcc")
#'
#' compare_results(res_methylcc, res_epidish)
#'
compare_results <- function(res1, res2) {
  if (is.null(res1) | is.null(res2)) {
    message("one of the inputs is NULL.")
    return()
  }
  if (!is.matrix(res1) | !is.matrix(res2)) {
    message("this function is for visualizing two of the this packages result matrices.")
    return()
  }
  
  df1 <- rename_cell_types(res1)
  df2 <- rename_cell_types(res2)
  
  df <- merge(df1, df2, by=0, all = T)
  df$Row.names <- NULL
  df <- df[, grepl(".", names(df), fixed = T)]
  CTs <- c()
  for (name in names(df)) {
    CT <- substr(name, 0, nchar(name)-2)
    if (!is.element(CT, CTs)) {
      CTs <- append(CTs, CT)
    }
  }
  df_final <- data.frame()
  for (CT in CTs) {
    col1 <- paste(CT, ".x", sep = "")
    col2 <- paste(CT, ".y", sep = "")
    df_temp <- df[, c(col1, col2)]
    names(df_temp)<- c("res1", "res2")
    df_temp <- cbind(df_temp, CT = c(CT))
    df_final <- rbind(df_final, df_temp)
  }
  ggplot2::ggplot(df_final, ggplot2::aes(x=res1, y=res2, color=CT)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline() +
    ggplot2::facet_wrap(~CT) +
    ggplot2::theme(legend.position = "none")
}
