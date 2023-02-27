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
    result$flowsortedbloodepic$B <- result$flowsortedbloodepic$Bcell
    result$flowsortedbloodepic$Bcell <- NULL
    result$flowsortedbloodepic$Neutro <- result$flowsortedbloodepic$Neu
    result$flowsortedbloodepic$Neu <- NULL
    df <- bind_rows(df, cbind(result$flowsortedbloodepic, method = c("FlowSortedBloodEPIC")))
  }
  if (!is.null(result$methylcc)) {
    result$methylcc$B <- result$methylcc$Bcell
    result$methylcc$Bcell <- NULL
    df <- bind_rows(df, cbind(result$methylcc, method = c("MethylCC")))
  }
  ggplot2::ggplot(data = reshape::melt(df), ggplot2::aes(x=variable,y=value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=variable)) +
    ggplot2::facet_wrap(~method) + labs(x = "Cell Type", y = "proportions") +
    ggplot2::theme(legend.position = "none")
}


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
  ggplot(df_final, aes(x=res1, y=res2, color=CT)) +
    geom_point() +
    geom_abline() +
    facet_wrap(~CT) +
    ggplot2::theme(legend.position = "none")
}
