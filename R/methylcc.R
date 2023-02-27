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
    return(NULL)
  }
  library(FlowSorted.Blood.450k)
  if (is(obj, "MethylSet")){
    obj <- minfi::mapToGenome(obj)
  }
  if (is(obj, "GenomicMethylSet")){
    set.seed(seed)
    return(as.data.frame(
      methylCC::cell_counts(methylCC::estimatecc(object = obj))))
  }
  if (!is(obj, "RGChannelSet")){
    message("input needs to be of class 'MethylSet', 'GenomicMethylSet' or 'RGChannelSet' to run methylCC deconvolution.")
    return(NULL)
  }
  set.seed(seed)
  return(as.data.frame(
    methylCC::cell_counts(methylCC::estimatecc(object = obj))))
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
    return(NULL)
  }
  if (all(is.na(unmeth))) {
    return(NULL)
  }
  methyl_set <- minfi::MethylSet(Meth = meth, Unmeth = unmeth)
  return(run_methylcc(methyl_set, seed = seed))
}
