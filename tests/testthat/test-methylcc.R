test_that("methylCC works", {
  requireNamespace("FlowSorted.Blood.450k", quietly = TRUE)
  requireNamespace("Biobase", quietly = TRUE)

  fsbe <- FlowSorted.Blood.450k::FlowSorted.Blood.450k
  rgset <- fsbe[,Biobase::pData(fsbe)$CellTypeLong %in% "Whole blood"]
  expect_equal(run_methylcc(rgset, 1), est)
})
