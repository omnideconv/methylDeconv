#' run EpiDISH
#'
#' @param value_matrix matrix containing methylation profiles
#'
#' @return named list with 3 Cell Type proportion matrices
#' @export
#'
#' @examples
run_epidish <- function(meth, unmeth, mode=c('RPC', 'CBS', 'CP'), 
                        reference=c('blood','breast','epithelial'), 
                        maxit = 50, nu.v = c(0.25, 0.5, 0.7), 
                        constraint = c("inequality", "equality")){
  require(EpiDISH)
  
  check_input(meth, unmeth)
  beta_matrix <- create_beta(meth, unmeth)
  
  if (length(mode) > 1) {
    mode <- mode[1]
    message(paste0(mode, " was chosen because multiple values were supplied for \"mode\""))
  }
  if (length(reference) > 1) {
    reference <- reference[1]
    message(paste0(reference, " was chosen because multiple values were supplied for \"reference\""))
  }
  message(paste0("Starting EpiDISH deconvolution with mode ", mode, " ..."))

  result_epidish <- switch (reference,
    'blood' = EpiDISH::epidish(beta.m = beta_matrix,
                               ref.m = EpiDISH::centDHSbloodDMC.m,
                               method = mode,
                               maxit = maxit,
                               nu.v = nu.v,
                               constraint = constraint),
    'breast' = EpiDISH::epidish(beta.m = beta_matrix,
                                ref.m = EpiDISH::centEpiFibFatIC.m,
                                method = mode,
                                maxit = maxit,
                                nu.v = nu.v,
                                constraint = constraint),
    'epithelial' = EpiDISH::epidish(beta.m = beta_matrix,
                                    ref.m = EpiDISH::centEpiFibFatIC.m,
                                    method = mode,
                                    maxit = maxit,
                                    nu.v = nu.v,
                                    constraint = constraint)
  )

  return(result_epidish)
}


