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
