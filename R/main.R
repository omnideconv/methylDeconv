#' Main function
#'
#' @return A String
#' @export
#'
#' @examples
run_deconvolutions <- function(obj=NA, meth=NA, unmeth=NA, meta=NA, condition=NA, seed=1){
  result_epidish <- NA
  result_tca <- NA
  result_flowsortedbloodepic <- NA
  result_methylcc <- NA
  if (missing(obj)) {
    if (missing(meth) | missing(unmeth)) {
      message("if no object is given, 'meth' and 'unmeth' need to be provided.")
      return(NA)
    }
    if (!is.matrix(meth)) {
      message("'meth' needs to be of class 'matix'.")
      return(NA)
    }
    if (!is.matrix(unmeth)) {
      message("'unmeth' needs to be of class 'matix'.")
      return(NA)
    }
    result_flowsortedbloodepic <- run_flowsortedbloodepic_raw(meth, unmeth)
    result_methylcc <- run_methylcc_raw(meth, unmeth)
  } else if (is.matrix(obj)) {
    result_epidish <- run_epidish(obj)
    if (!missing(meta) & !missing(condition)) {
      if (!is.data.frame(meta) | !is.character(condition)) {
        message("'meta' needs to be of class 'data.frame' and 'condition' of class 'character'.")
      } else {
        result_tca <- run_tca(obj, meta, condition, result_epidish$rpc)
      }
    }
    if (!missing(meth) & !missing(unmeth)) {
      if (!is.matrix(meth)) {
        message("'meth' needs to be of class 'matix'.")
        return(NA)
      }
      if (!is.matrix(unmeth)) {
        message("'unmeth' needs to be of class 'matix'.")
        return(NA)
      }
      result_flowsortedbloodepic <- run_flowsortedbloodepic_raw(meth, unmeth)
      result_methylcc <- run_methylcc_raw(meth, unmeth, seed)
    }
  } else if (class(obj) == "MethylSet") {
    result_flowsortedbloodepic <- run_flowsortedbloodepic(obj)
    result_methylcc <- run_methylcc(obj, seed)
    if (!missing(meth) & !missing(unmeth)) {
      message("'meth' and 'unmeth' will be ignored, since 'obj' is of class 'MethylSet'.")
    }
  } else if (class(obj) == "GenomicMethylSet") {
    result_flowsortedbloodepic <- run_flowsortedbloodepic(obj)
    result_methylcc <- run_methylcc(obj, seed)
    if (!missing(meth) & !missing(unmeth)) {
      message("'meth' and 'unmeth' will be ignored, since 'obj' is of class 'GenomicMethylSet'.")
    }
  } else if (class(obj) == "RGChannelSet") {
    result_flowsortedbloodepic <- run_flowsortedbloodepic(obj, method = "preprocessNoob")
    result_methylcc <- run_methylcc(obj, seed)
    if (!missing(meth) & !missing(unmeth)) {
      message("'meth' and 'unmeth' will be ignored, since 'obj' is of class 'RGChannelSet'.")
    }
  } else {
    message("'obj' needs to be of class 'matrix', 'MethylSet', 'GenomicMethylSet' or 'RGChannelSet'.")
    return(NA)
  }
  return(list(epidish=result_epidish,
              tca=result_tca,
              flowsortedbloodepic=result_flowsortedbloodepic,
              methylcc=result_methylcc))
}
