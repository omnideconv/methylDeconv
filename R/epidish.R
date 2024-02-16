#' run EpiDISH
#'
#' @param meth methylated data matrix
#' @param unmeth unmethylated data matrix
#' @param mode Choice of a reference-based method ('RPC','CBS','CP')
#' @param reference A matrix of reference 'centroids', i.e. representative molecular profiles, 
#' for a number of cell subtypes. rows label molecular features (e.g. CpGs,...) 
#' and columns label the cell-type. IDs need to be provided as rownames and 
#' colnames, respectively. Missing value is not allowed, and all values in 
#' this matrix should be positive or zero. For DNAm data, values should be 
#' beta-values.
#' @param maxit Only used in RPC mode, the limit of the number of IWLS iterations
#' @param nu.v Only used in CBS mode. It is a vector of several candidate nu values. nu is 
#' parameter needed for nu-classification, nu-regression, and 
#' one-classification in svm. The best estimation results among all candidate nu 
#' will be automatically returned.
#' @param constraint Only used in CP mode, you can choose either of 'inequality' or 'equality' 
#' normalization constraint. The default is 'inequality' (i.e sum of weights 
#' adds to a number less or equal than 1), which was implemented in 
#' Houseman et al (2012).
#'
#' @return CP-mode
#' A list with the following entries: estF: a matrix of the estimated fractions; 
#' ref: the reference centroid matrix used; dataREF: the subset of the input 
#' data matrix with only the probes defined in the reference matrix.
#' 
#' @return CBS-mode
#' A list with the following entries: estF: a matrix of the estimated fractions; 
#' nu: a vector of 'best' nu-parameter for each sample; 
#' ref: the reference centroid matrix used;
#' dataREF: the subset of the input data matrix with only the probes defined in the 
#' reference matrix.
#' 
#' @return RPC-mode
#' A list with the following entries: estF: a matrix of the estimated fractions;
#'  ref: the reference centroid matrix used; 
#' dataREF: the subset of the input data matrix with only the probes defined in the 
#' reference matrix.
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
    message(paste0(mode, " was chosen as default for \"mode\""))
  }
  if (length(reference) > 1) {
    reference <- reference[1]
    message(paste0(reference, " was chosen as default for \"reference\""))
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
                                    ref.m = EpiDISH::centEpiFibIC.m,
                                    method = mode,
                                    maxit = maxit,
                                    nu.v = nu.v,
                                    constraint = constraint)
  )

  return(result_epidish)
}


