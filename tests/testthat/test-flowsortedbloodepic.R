test_that("FlowSortedBloodEPIC works", {
  FlowSorted.Blood.EPIC <- FlowSorted.Blood.EPIC::libraryDataGet('FlowSorted.Blood.EPIC')
  RGsetTargets <- FlowSorted.Blood.EPIC[,FlowSorted.Blood.EPIC$CellType == "MIX"]
  sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
                                     seq_len(dim(RGsetTargets)[2]),
                                     sep = "_")
  expect_equal(run_flowsortedbloodepic(RGsetTargets, method="preprocessNoob"), propepic)
})
