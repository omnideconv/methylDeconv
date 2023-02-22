#' run FlowSortedBloodEPIC deconvolution
#'
#' @param methyl_set minfi::MethylSet
#'
#' @return Cell Type proportion matrix
#' @export
#'
#' @examples
run_flowsortedbloodepic <- function(obj, method="preprocessQuantile"){
  if (length(obj) == 1) {
    return(NA)
  }
  return(FlowSorted.Blood.EPIC::estimateCellCounts2(obj,
                                                    compositeCellType = "Blood",
                                                    processMethod = method,
                                                    probeSelect = "IDOL",
                                                    cellTypes = c("CD8T",
                                                                  "CD4T",
                                                                  "NK",
                                                                  "Bcell",
                                                                  "Mono",
                                                                  "Neu"))$prop)
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
  if (all(is.na(meth))) {
    return(NA)
  }
  if (all(is.na(unmeth))) {
    return(NA)
  }
  methyl_set <- minfi::MethylSet(Meth = meth, Unmeth = unmeth)
  return(run_flowsortedbloodepic(methyl_set))
}
