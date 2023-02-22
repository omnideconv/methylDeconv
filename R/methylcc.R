#' run methylCC deconvolution
#'
#' @param genomic_methyl_set minfi::MethylSet or minfi::GenomicMethylSet
#' @param seed numeric
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_methylcc <- function(obj, seed=1){
  if (length(obj) == 1) {
    return(NA)
  }
  library(FlowSorted.Blood.450k)
  if (class(obj) == "MethylSet"){
    obj <- minfi::mapToGenome(obj)
  }
  if (class(obj) == "GenomicMethylSet"){
    set.seed(seed)
    return(methylCC::cell_counts(methylCC::estimatecc(object = obj)))
  }
  if (class(obj) != "RGChannelSet"){
    message("input needs to be of class 'MethylSet', 'GenomicMethylSet' or 'RGChannelSet' to run methylCC deconvolution.")
    return(NA)
  }
  set.seed(seed)
  return(methylCC::cell_counts(methylCC::estimatecc(object = obj)))
}

#' run methylCC deconvolution using raw files
#'
#' @param meth raw methylation data (meth)
#' @param unmeth raw methylation data (unmeth)
#' @param seed numeric
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_methylcc_raw <- function(meth, unmeth, seed = 1){
  if (all(is.na(meth))) {
    return(NA)
  }
  if (all(is.na(unmeth))) {
    return(NA)
  }
  methyl_set <- minfi::MethylSet(Meth = meth, Unmeth = unmeth)
  return(run_methylcc(methyl_set, seed = seed))
}
