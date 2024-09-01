
#' List of supported deconvolution methods
#'
#' The methods currently supported are
#' `EpiDISH`, `FlowSorted`, `MethylCC`, `MethylResolver`, `MethAtlas`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods <- c(
  "EpiDISH" = "epidish", "FlowSorted" = "flowsorted", "MethylCC" = "methylcc", 
  "MethylResolver" = "methylresolver", "MethAtlas" = "meth_atlas"
)


#' Deconvolution
#'
#' @param methyl_set A minfi MethylSet
#' @param method A string specifying the method. Supported methods are 'epidish', 'flowsorted', 'methylcc', 'methylresolver'
#' @param normalize_results   Whether the deconvolution results should be normalized.
#'   Negative values will be put to 0, and the estimates will be normalized to sum to 1.
#'   Defaults to FALSE.
#' @param ... Additional parameters, passed to the algorithm used. See individual method documentations for details.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are
#' individuals, columns are cell types.
#' @export
#'
#' @examples 
#' 
#' ex_data <- minfiData::MsetEx
#' 
#' result <- deconvolute(ex_data, method='epidish')
deconvolute <- function(methyl_set, method=deconvolution_methods, normalize_results = FALSE, ...){
  
  if (length(method) > 1) {
    stop(
      "Please only specify one method and not ", length(method), ": ",
      paste(method, collapse = ", ")
    )
  }
  
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  
  method <- tolower(method)
  
  
  result <- switch (method,
    epidish = run_epidish(methyl_set, ...)$estF,
    flowsorted = run_flowsortedblood(methyl_set, ...)$prop,
    methylcc = as.matrix(run_methylcc(methyl_set, ...)),
    methylresolver = as.matrix(run_methylresolver(methyl_set, ...)$result_fractions),
    meth_atlas = run_meth_atlas(methyl_set, ...)
  )
  
  if(!is.null(result)){
    # Normalize the results to sum up to 1
    if (normalize_results) {
      deconv <- normalize_deconv_results(result)
    }
    # Alphabetical order of celltypes
    result <- result[, order(colnames(result)), drop = FALSE]
  }
  
  return(result)
}

#' Run all available deconvolution methods
#'
#' @param methyl_set A minfi MethylSet
#' @param array type of methylation array that was used. possible options are '450k' and 'EPIC'
#'
#' @return dataframe with results of all methods
#' @export
#'
#' @examples
run_all_methods <- function(methyl_set, array = c('450k','EPIC')){
  
  res_epidish <- run_epidish(methyl_set)
  res_flowsorted <- run_flowsortedblood(methyl_set, array = array)
  res_methylcc <- run_methylcc(methyl_set, array = array)
  res_methylresolver <- run_methylresolver(methyl_set)
  res_meth_atlas <- run_meth_atlas(methyl_set)
  
  results <- list(res_epidish$estF, res_flowsorted$prop, res_methylcc, res_methylresolver$result_fractions, res_meth_atlas)
  names(results) <- c('EpiDISH', "FlowSortedBlood", "MethylCC", "MethylResolver", "MethAtlas")
  
  tmp <- lapply(1:4, function(i){
    result_i <- results[[i]]
    print(result_i)
    result_i <- data.frame(result_i[, order(colnames(result_i)), drop = FALSE], check.names = F)
    
    result_i <- rename_cell_types(result_i)
    result_i <- result_i[,colnames(result_i) != "other"]
    print(result_i)
    
    result_i <- tibble::rownames_to_column(result_i, "sample")
    result_i[['method']] <- names(results)[i]
    
    print(result_i)
    result_i
  })
  
  combined_result <- data.table::rbindlist(tmp, fill = TRUE)
  
  return(combined_result)
}
