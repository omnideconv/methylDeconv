% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylcc.R
\name{run_methylcc}
\alias{run_methylcc}
\title{run methylCC deconvolution}
\usage{
run_methylcc(
  methyl_set,
  array = c("450k", "EPIC"),
  find_dmrs_object = NULL,
  verbose = TRUE,
  epsilon = 0.01,
  max_iter = 100,
  take_intersection = FALSE,
  include_cpgs = FALSE,
  include_dmrs = TRUE,
  init_param_method = "random",
  a0init = NULL,
  a1init = NULL,
  sig0init = NULL,
  sig1init = NULL,
  tauinit = NULL,
  demo = FALSE,
  seed = 1
)
}
\arguments{
\item{methyl_set}{A minfi MethylSet}

\item{array}{type of methylation array that was used. possible options are '450k' and 'EPIC'}

\item{find_dmrs_object}{If the user would like to supply different differentially methylated regions, they can use the output from the find_dmrs function to supply different regions to estimatecc.}

\item{verbose}{TRUE/FALSE argument specifying if verbose messages should be returned or not. Default is TRUE.}

\item{epsilon}{Threshold for EM algorithm to check for convergence. Default is 0.01.}

\item{max_iter}{Maximum number of iterations for EM algorithm. Default is 100 iterations.}

\item{take_intersection}{TRUE/FALSE asking if only the CpGs included in object should be used to find DMRs. Default is FALSE.}

\item{include_cpgs}{TRUE/FALSE. Should individual CpGs be returned. Default is FALSE.}

\item{include_dmrs}{TRUE/FALSE. Should differentially methylated regions be returned. Default is TRUE.}

\item{init_param_method}{method to initialize parameter estimates. Choose between "random" (randomly sample) or "known_regions" (uses unmethyalted and methylated regions that were identified based on Reinus et al. (2012) cell sorted data.). Defaults to "random".}

\item{a0init}{Default NULL. Initial mean methylation level in unmethylated regions}

\item{a1init}{Default NULL. Initial mean methylation level in methylated regions}

\item{sig0init}{Default NULL. Initial var methylation level in unmethylated regions}

\item{sig1init}{Default NULL. Initial var methylation level in methylated regions}

\item{tauinit}{Default NULL. Initial var for measurement error}

\item{demo}{TRUE/FALSE. Should the function be used in demo mode to shorten examples in package. Defaults to FALSE.}

\item{seed}{fixed seed to account for RNG influences}
}
\value{
A object of the class estimatecc that contains information about the cell composition estimation (in the summary slot) and the cell composition estimates themselves (in the cell_counts slot).
}
\description{
run methylCC deconvolution
}
