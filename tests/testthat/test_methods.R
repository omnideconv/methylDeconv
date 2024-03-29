library(methylDeconv)
library(minfiData)
library(minfi)

ex_data <- minfiData::MsetEx
meth <- minfi::getMeth(ex_data)
unmeth <- minfi::getUnmeth(ex_data)

# write.csv()

test_that("EpiDISH works", {
  epidish_res <- methylDeconv::run_epidish(meth = meth, 
                                           unmeth = unmeth)$estF
  check_result <- as.matrix(read.csv("test_results/epidish.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = epidish_res,
    expected = check_result, tolerance = 1e-3
  )
})


test_that("methylCC works", {
  methylcc_res <- as.matrix(methylDeconv::run_methylcc(meth = meth, 
                                                       unmeth = unmeth))
  check_result <- as.matrix(read.csv("test_results/methylcc.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = methylcc_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("FlowSorted works", {
  flowSorted_res <- methylDeconv::run_flowsortedblood(meth = meth, 
                                                      unmeth = unmeth)$prop
  check_result <- as.matrix(read.csv("test_results/flowsorted.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = flowSorted_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("MethylResolver works", {
  methylResolver_res <- methylDeconv::run_methylresolver(meth = meth, 
                                                         unmeth = unmeth)$result_fractions
  check_result <- as.matrix(read.csv("test_results/methylresolver.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = methylResolver_res,
    expected = check_result, tolerance = 1e-3
  )
})
