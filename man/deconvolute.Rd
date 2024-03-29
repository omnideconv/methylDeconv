% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{deconvolute}
\alias{deconvolute}
\title{Deconvolution}
\usage{
deconvolute(
  meth,
  unmeth,
  method = deconvolution_methods,
  normalize_results = FALSE,
  ...
)
}
\arguments{
\item{meth}{methylated data matrix}

\item{unmeth}{unmethylated data matrix}

\item{method}{A string specifying the method. Supported methods are 'epidish', 'flowsorted', 'methylcc', 'methylresolver'}

\item{normalize_results}{Whether the deconvolution results should be normalized.
Negative values will be put to 0, and the estimates will be normalized to sum to 1.
Defaults to FALSE.}

\item{...}{Additional parameters, passed to the algorithm used. See individual method documentations for details.}
}
\value{
A matrix with the probabilities of each cell-type for each individual. Rows are
individuals, columns are cell types.
}
\description{
Deconvolution
}
\examples{

ex_data <- minfiData::MsetEx
meth <- minfi::getMeth(ex_data)
unmeth <- minfi::getUnmeth(ex_data)

result <- deconvolute(meth, unmeth, method='epidish')
}
