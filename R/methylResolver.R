#' Run MethylResolver
#'
#' @param methyl_set A minfi MethylSet
#' @param doPar Whether to use parallel processing to speed up the deconvolution computation if many samples are present. Default is FALSE. 
#' @param numCores Number of cores used for parallel processing to speed up the deconvolution of many samples. Requires doPar = TRUE. Default is 1. numCores = "auto" is max number of cores available minus one. 
#' @param alpha Set the alpha parameter for LTS deconvolution. This is the fraction of optimal CpGs from the signature matrix which are used for deconvolution. Must be between 0 and 1. Users can specify a vector or a single number. If a vector is specified, a grid search of the values is conducted and the alpha value that results in the lowest RMSE between the original and reconstructed mixture is selected. Default is seq(0.5,0.9,by = 0.05). 
#' @param absolute Whether to compute tumor purity and absolute cell type fractions. Default is TRUE. 
#' @param purityModel Random Forest model to predict mixture purity (unknown content) which allows the calculation of absolute cell type fractions. Required if absolute is TRUE. Default is our RF model trained on the consensus purity estimate (CPE) using TCGA data.  
#' @param seed fixed seed to account for RNG influences
#'
#' @return 
#' @export
#'
#' @examples
run_methylresolver <- function(methyl_set, doPar = F, numCores = 1, alpha = seq(0.5,0.9,by = 0.05),
                               absolute = TRUE, purityModel = MethylResolver::RFmodel, seed = 1){
  
  set.seed(seed)
  
  check_input_mset(methyl_set) 
  beta_matrix <- minfi::getBeta(methyl_set)
  
  if (length(alpha) > 1){
    warning("MethylResolver may fail if multiple alpha values are provided. If this occurs, specify a single alpha value between 0.5 and 1.",
            immediate. = TRUE)
  }
  
  result_methylresolver <- MethylResolver::MethylResolver(methylMix = beta_matrix, 
                                                          methylSig = MethylResolver::MethylSig, 
                                                          betaPrime = FALSE,
                                                          doPar = doPar, 
                                                          numCores = numCores, 
                                                          alpha = alpha, 
                                                          absolute = absolute, 
                                                          purityModel = purityModel)
  
  result_metrics <- result_methylresolver[,1:4]
  result_fractions <- result_methylresolver[,5:15]
  result_absolute <- result_methylresolver[,16:26]
  result_purity <- result_methylresolver[,27]
  
  return(list(result_metrics=result_metrics,
              result_fractions=result_fractions,
              result_absolute=result_absolute,
              result_purity=result_purity))
}