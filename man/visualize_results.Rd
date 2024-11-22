% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{visualize_results}
\alias{visualize_results}
\title{Function to plot the results from \code{deconvolute()} or \code{run_all_methods()} as boxplots.}
\usage{
visualize_results(result)
}
\arguments{
\item{result}{result from \code{deconvolute()} or \code{run_all_methods()}}
}
\description{
Function to plot the results from \code{deconvolute()} or \code{run_all_methods()} as boxplots.
}
\examples{

mset <- minfiData::MsetEx

# one result
res_epidish <- deconvolute(mSet, "epidish")
visualize_results(res_epidish)

# all results
results <- run_all_methods(mSet)
visualize_results(res_epidish)

# combine mutliple results into a list
res_methylcc <- deconvolute(mset, "methylcc")
visualize_results(list(epidish = res_epidish, 
                       methycc = res_methylcc))

}
