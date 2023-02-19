#' run EpiDISH
#'
#' @param value_matrix matrix containing methylation profiles
#'
#' @return named list with 3 Cell Type proportion matrices
#' @export
#'
#' @examples
run_epidish <- function(value_matrix){
  message("running EpiDISH deconvolutions:")
  rpc_result <- run_epidish_rpc(value_matrix)
  cbs_result <- run_epidish_cbs(value_matrix)
  cp_result <- run_epidish_cp(value_matrix)
  return(list(rpc=rpc_result, cbs=cbs_result, cp=cp_result))
}

#' Title
#'
#' @param value_matrix matrix containing methylation profiles
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_epidish_rpc <- function(value_matrix){
  message("1. using 'Robust Partial Correlations'.")
  return(EpiDISH::epidish(beta.m = value_matrix,
                 ref.m = EpiDISH::centDHSbloodDMC.m,
                 method = 'RPC')
         $estF)
}

#' Title
#'
#' @param value_matrix matrix containing methylation profiles
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_epidish_cbs <- function(value_matrix){
  message("2. using 'CIBERSORT'.")
  return(EpiDISH::epidish(beta.m = value_matrix,
                 ref.m = EpiDISH::centDHSbloodDMC.m,
                 method = 'CBS')
         $estF)
}

#' Title
#'
#' @param value_matrix matrix containing methylation profiles
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_epidish_cp <- function(value_matrix){
  message("3. using 'Constrained Projection '.")
  return(EpiDISH::epidish(beta.m = value_matrix,
                 ref.m = EpiDISH::centDHSbloodDMC.m,
                 method = 'CP')
         $estF)
}
