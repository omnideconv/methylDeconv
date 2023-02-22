test_that("main works", {
  # no input
  expect_no_error(run_deconvolutions())
  # only matrix
  expect_no_error(run_deconvolutions(obj=matrix()))
  # matrix and meta + condition
  expect_no_error(run_deconvolutions(obj=matrix(), meta=data.frame(), condition =""))
  # MethylSet, GenomicMethylSet, RGChannelSet
  expect_no_error(run_deconvolutions(obj=minfi::MethylSet(Meth=matrix(), Unmeth=matrix())))
  # wrong obj class
  expect_no_error(run_deconvolutions(obj="wrong type"))
  # meta without condition
  expect_no_error(run_deconvolutions(obj=matrix(), meta=data.frame()))
  # condition without meta
  expect_no_error(run_deconvolutions(obj=matrix(), condition="condition"))
  # meth without unmeth
  expect_no_error(run_deconvolutions(meth=matrix()))
  # unmeth without meth
  expect_no_error(run_deconvolutions(obj=matrix(), unmeth=matrix()))
})
