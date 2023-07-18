test_that("TCA works", {
  test_result <- run_tca(hannum$X, hannum$cov, c("gender","age"), hannum$W)
  expect_equal(test_result, tca.mdl.hannum[["W"]])
})
