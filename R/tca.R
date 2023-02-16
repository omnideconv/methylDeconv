#' run TCA
#'
#' @param value_matrix matrix containing methylation profiles
#' @param meta_data matrix containing information about samples in 'value_matrix'
#' @param condition collumn name of 'meta_data' either complex or numeric
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_tca <- function(value_matrix, meta_data, condition){
  message("running TCA deconvolution.")
  meta_data <-  meta_data[,condition, drop = FALSE]
  if (is.na(value_matrix)){
    message("ommited rows containing missing values for TCA analysis")
    value_matrix <- na.omit(value_matrix)
  }
  if (sum(rowVars(value_matrix) < 1e-08) != 0){
    message("removing rows with variance less than 1e-08 for TCA analysis")
    value_matrix <- value_matrix[rowVars(value_matrix) >= 1e-08,]
  }
  tca_result <- TCA::tca(X = value_matrix,
                      W = TCA::BloodFrac.RPC.m,
                      C1 = meta_data,
                      )
  return(tca_result[["W"]])
}
