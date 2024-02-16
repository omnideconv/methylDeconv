#' Check input matrices
#'
#' @param meth 
#' @param unmeth 
#'
#'
#' @examples
check_input <- function(meth, unmeth){
  if(!all(dim(meth) == dim(unmeth))){
    stop('Meth and Unmeth matrices have different dimensions, stopping.')
  }
  if(all(is.na(meth))){
    stop('All entries in Meth matrix are NA, stopping.')
  }
  if(all(is.na(unmeth))){
    stop('All entries in Unmeth matrix are NA, stopping.')
  }
}


#' Create beta matrix from meth and unmeth matrices
#'
#' @param meth 
#' @param unmeth 
#'
#' @return beta matrix
#'
#' @examples
create_beta <- function(meth, unmeth){
  methyl_set <- minfi::MethylSet(Meth = meth, Unmeth = unmeth)
  ratio_set <- minfi::ratioConvert(methyl_set)
  
  return(minfi::getBeta(ratio_set))
}


#' Normalize deconvolution result
#'
#' @param deconv_result The original deconvolution result
#'
#' @return A matrix with the rowsums of one and no negative values
#' @export
normalize_deconv_results <- function(deconv_result) {
  celltypes <- colnames(deconv_result)
  deconv_result[deconv_result < 0] <- 0
  deconv_result <- t(apply(deconv_result, 1, function(row) row / sum(row)))
  # Apply returns a vector when only supplied with one celltype. To counter it and return a matrix
  # and not a vector, this operation is needed
  if (length(celltypes) == 1) {
    deconv_result <- t(deconv_result)
    colnames(deconv_result) <- celltypes
  }
  return(deconv_result)
}