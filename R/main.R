
#' List of supported deconvolution methods
#'
#' The methods currently supported are
#' `EpiDISH`, `Houseman`, `MethylCC`, `MethylResolver`, `MethAtlas`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods <- c(
  "EpiDISH" = "epidish", "Houseman" = "houseman", "MethylCC" = "methylcc", 
  "MethylResolver" = "methylresolver", "MethAtlas" = "methatlas"
)


#' Deconvolution with methyldeconv
#'
#' @param methyl_set A minfi MethylSet
#' @param method A string specifying the method. Supported methods are 'epidish', 'houseman', 'methylcc', 'methylresolver', 'methatlas'
#' @param scale_results   Whether the deconvolution results should be rescaled.
#'   Negative values will be set to 0, and the estimates will be normalized to sum to 1 per sample.
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
deconvolute <- function(methyl_set, method=deconvolution_methods, scale_results = FALSE, ...){
  
  options(matrixStats.useNames.NA = "deprecated")
  
  if (length(method) > 1) {
    stop(
      "Please only specify one method and not ", length(method), ": ",
      paste(method, collapse = ", ")
    )
  }
  
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  
  check_input_mset(methyl_set)
  beta_matrix <- create_beta(methyl_set)
  
  method <- tolower(method)
  
  result <- switch (method,
    epidish = run_epidish(beta_matrix, ...)$estF,
    houseman = run_houseman(methyl_set, ...)$prop,
    methylcc = as.matrix(run_methylcc(methyl_set, ...)),
    methylresolver = as.matrix(run_methylresolver(beta_matrix, ...)$result_fractions),
    methatlas = run_methatlas(beta_matrix, ...)
  )
  
  if(!is.null(result)){
    # Scale the results to sum up to 1
    if (scale_results) {
      deconv <- normalize_deconv_results(result)
    }
    # Alphabetical order of celltypes
    result <- result[, order(colnames(result)), drop = FALSE]
  }
  
  return(result)
}

#' Run selected set of methods and aggregate results
#'
#' @param methyl_set A minfi MethylSet
#' @param array type of methylation array that was used. possible options are '450k' and 'EPIC'
#' @param methods list of methods (>1) that will be applied to the methyl set
#' @param aggregation_method how should deconvolution estimates be combined (default is mean)
#' @param scale_results   Whether the deconvolution results should be rescaled.
#'   Negative values will be set to 0, and the estimates will be normalized to sum to 1 per sample.
#'   Defaults to FALSE.
#'   
#' @return dataframe with results of all selected methods as well as the combined estimates
#' @export
#'
deconvolute_combined <- function(methyl_set, array = c('450k','EPIC'), methods, aggregation_method = 'mean', scale_results = FALSE){
  
  if(any(!method %in% deconvolution_methods)){
    stop(paste0('At least one of your selected methods is not supported by methyldeconv. Please check your spelling, supported methods are: ',
                'epidish, houseman, methylcc, methylresolver, methatlas'))
  }
  
  result <- lapply(methods, function(m){
    deconvolute(methyl_set = methyl_set, method = m, scale_results = scale_results)
  })
  
  # TODO complete this
  
  
  methyldeconv_df <- methyldeconv_df |> 
    group_by(sample, celltype_clean) |> 
    dplyr::summarize(mean_value = mean(value), .groups = 'drop') |>
    group_by(sample) |> 
    do(mutate(., value = mean_value / sum(mean_value))) |>
    ungroup() |>
    mutate(method = 'combined') |> 
    bind_rows(methyldeconv_df) 


  # for methylresolver, combine Tnaive and Tmem to CD4+ cells
  res_methylresolver <- as.matrix(res_methylresolver$result_fractions)
  res_methylresolver <- cbind(res_methylresolver, 
                              "T cell CD4+" = res_methylresolver[, "Tmem"] + res_methylresolver[, "Tnaive"])
  
  results <- list(res_epidish$estF, res_flowsorted$prop, res_methylcc, res_methylresolver, res_meth_atlas)
  names(results) <- c('EpiDISH', "Houseman", "MethylCC", "MethylResolver", "MethAtlas")
  
  tmp <- lapply(1:5, function(i){
    result_i <- results[[i]]
    result_i <- data.frame(result_i[, order(colnames(result_i)), drop = FALSE], check.names = F)
    
    result_i <- rename_cell_types(result_i)
    result_i <- result_i[,colnames(result_i) != "other"]
    
    result_i <- tibble::rownames_to_column(result_i, "sample")
    result_i[['method']] <- names(results)[i]
    
    result_i
  })
  
  combined_result <- data.table::rbindlist(tmp, fill = TRUE)
  
  return(combined_result)
}
