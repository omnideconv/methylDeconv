test_that("EpiDISH works", {
  test_output <- run_epidish(EpiDISH::DummyBeta.m)
  expect_equal(test_output$rpc, Bloodfrac.RPC.m)
  expect_equal(test_output$cbs, Bloodfrac.CBS.m)
  expect_equal(test_output$cp, Bloodfrac.CP.m)
})
