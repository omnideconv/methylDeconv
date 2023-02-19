#' run FlowSortedBloodEPIC deconvolution
#'
#' @param methyl_set minfi::MethylSet
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_flowsortedbloodepic <- function(methyl_set, method="preprocessQuantile"){
  return(FlowSorted.Blood.EPIC::estimateCellCounts2(methyl_set,
                                                    compositeCellType = "Blood",
                                                    processMethod = method,
                                                    probeSelect = "IDOL",
                                                    cellTypes = c("CD8T",
                                                                  "CD4T",
                                                                  "NK",
                                                                  "Bcell",
                                                                  "Mono",
                                                                  "Neu")))$prop
}

#' run FlowSortedBloodEPIC deconvolution using raw files
#'
#' @param meth raw methylation data (meth)
#' @param unmeth raw methylation data (unmeth)
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_flowsortedbloodepic_raw <- function(meth, unmeth){
  methyl_set <- minfi::MethylSet(Meth = meth, Unmeth = unmeth)
  return(run_flowsortedbloodepic(methyl_set))
}
