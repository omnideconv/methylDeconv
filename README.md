# methylDeconv

<!-- badges: start -->

[![R-CMD-check](https://github.com/omnideconv/methylDeconv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/omnideconv/methylDeconv/actions/workflows/R-CMD-check.yml) [![Codecov test coverage](https://codecov.io/gh/omnideconv/methylDeconv/branch/main/graph/badge.svg)](https://app.codecov.io/gh/omnideconv/methylDeconv?branch=main)

<!-- badges: end -->

Ever wanted to apply cell-type deconvolution on your DNA methylation data but could not decide which method to use? Here is **methylDeconv** to help your needs!

This package integrates unified access to five reference-based cell-type deconvolution methods that can directly be applied to Illumina array data (450k, EPIC arrays) or bisulfite sequencing data (RRBS, WGBS).

The included methods are:

| method                       | license | citation |
|------------------------------|---------|----------|
| EpiDISH                      |         |          |
| Houseman (Flow.Sorted.Blood) |         |          |
| MethAtlas                    |         |          |
| methylCC                     |         |          |
| methylResolver               |         |          |

## Installation

You can install methylDeconv from [GitHub](https://github.com/), we recommend to use the [pak](https://github.com/r-lib/pak) package manager:

``` r
# install the `pak` package manager
install.packages("pak")

pak::pkg_install("omnideconv/methylDeconv")
```

## Example

methylDeconv can either be applied directly to a methylSet from the minfi package, or you can apply each method separately on a beta matrix with Illumina CpG IDs.

Both cases will be demonstrated here using example data from minfi:

``` r
library(methylDeconv)
## basic example code
```
