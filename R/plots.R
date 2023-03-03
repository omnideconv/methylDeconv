#' Title
#'
#' @param result list of 4
#'
#' @return
#' @export
#'
#' @examples
visualize_results <- function(result){
  if (!is.list(result)) {
    message("this function is for visualizing the this packages result object.")
    return()
  }
  df <- data.frame()
  if (!is.null(result$epidish)) {
    if (!is.null(result$epidish$rpc)) {
      df <- rbind(df, cbind(result$epidish$rpc, method = c("EpiDISH-RPC")))
    }
    if (!is.null(result$epidish$cbs)) {
      df <- rbind(df, cbind(result$epidish$cbs, method = c("EpiDISH-CBS")))
    }
    if (!is.null(result$epidish$cp)) {
      df <- rbind(df, cbind(result$epidish$cp, method = c("EpiDISH-CP")))
    }
  }
  if (!is.null(result$tca)) {
    df <- rbind(df, cbind(result$tca, method = c("TCA")))
  }
  if (!is.null(result$flowsortedbloodepic)) {
    df <- dplyr::bind_rows(df, cbind(result$flowsortedbloodepic, method = c("FlowSortedBloodEPIC")))
  }
  if (!is.null(result$methylcc)) {
    df <- dplyr::bind_rows(df, cbind(result$methylcc, method = c("MethylCC")))
  }
  ggplot2::ggplot(data = reshape::melt(df), ggplot2::aes(x=variable,y=value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=variable)) +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(x = "Cell Type", y = "proportions") +
    ggplot2::theme(legend.position = "none")
}


#' Title
#'
#' @param result
#' @param CT
#'
#' @return
#' @export
#'
#' @examples
visualize_result_box <- function(result, CT='B') {
  if (!is.data.frame(result)) {
    message("this function is for visualizing one of the this packages result data frames.")
    return()
  }
  df <- as.data.frame(result[,CT])
  colnames(df) <- c("value")
  df <- cbind(df, CT = c(CT))
  ggplot2::ggplot(data = df, ggplot2::aes(x=CT,y=value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=CT)) +
    ggplot2::labs(y = "proportions") +
    ggplot2::theme(legend.position = "none")
}


#' Title
#'
#' @param result
#'
#' @return
#' @export
#'
#' @examples
visualize_result_bar <- function(result) {
  if (!is.data.frame(result)) {
    message("this function is for visualizing one of the this packages result data frames.")
    return()
  }
  ggplot2::ggplot(df, ggplot2::aes(x = sample, y = value, fill = variable)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::xlab("sample") +
    ggplot2::ylab("Cell Type probability distribution")
}


#' Title
#'
#' @param res1
#' @param res2
#'
#' @return
#' @export
#'
#' @examples
compare_results <- function(res1, res2) {
  if (is.null(res1) | is.null(res2)) {
    message("one of the inputs is NULL.")
    return()
  }
  df <- merge(res1, res2, by=0, all = T)
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
