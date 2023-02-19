#' run methylCC deconvolution
#'
#' @param genomic_methyl_set minfi::MethylSet or minfi::GenomicMethylSet
#' @param seed numeric
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_methylcc <- function(genomic_methyl_set, seed = 1){
  if (class(genomic_methyl_set) == "MethylSet"){
    genomic_methyl_set <- minfi::mapToGenome(genomic_methyl_set)
  }
  if (class(genomic_methyl_set) == "GenomicMethylSet"){
    set.seed(seed)
    return(methylCC::cell_counts(methylCC::estimatecc(object = genomic_methyl_set)))
  }
  if (class(genomic_methyl_set) != "RGChannelSet"){
    message("input needs to be of class 'MethylSet', 'GenomicMethylSet' or 'RGChannelSet' to run methylCC deconvolution.")
    return(NA)
  }
  set.seed(seed)
  est <- methylCC::estimatecc(object = genomic_methyl_set)
  return(methylCC::cell_counts(est))
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
  methyl_set <- minfi::MethylSet(Meth = meth, Unmeth = unmeth)
  return(run_methylcc(methyl_set, seed = seed))
}
