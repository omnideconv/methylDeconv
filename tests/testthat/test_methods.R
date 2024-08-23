library(methylDeconv)
library(minfiData)
library(minfi)

ex_data <- minfiData::MsetEx
meth <- minfi::getMeth(ex_data)
unmeth <- minfi::getUnmeth(ex_data)

# write.csv()

test_that("EpiDISH works", {
  epidish_res <- methylDeconv::run_epidish(methyl_set = ex_data)$estF
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
  methylcc_res <- as.matrix(methylDeconv::run_methylcc(methyl_set = ex_data))
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
  flowSorted_res <- methylDeconv::run_flowsortedblood(methyl_set = ex_data)$prop
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
  methylResolver_res <- as.matrix(methylDeconv::run_methylresolver(methyl_set = ex_data, alpha = 1)$result_fractions)
  check_result <- as.matrix(read.csv("test_results/methylresolver.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = methylResolver_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("MethAtlas works", {
  meth_atlas_res <- methylDeconv::run_meth_atlas(meth = meth, 
                                           unmeth = unmeth)
  check_result <- as.matrix(read.csv("test_results/meth_atlas.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = meth_atlas_res,
    expected = check_result, tolerance = 1e-3
  )
})
