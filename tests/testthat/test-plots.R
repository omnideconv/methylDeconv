test_that("plotting works", {
  # wrong input
  expect_no_error(visualize_results(0))
  # correct result
  expect_no_error(visualize_results(result))
  # wrong input
  expect_no_error(visualize_result_box(0))
  # correct result
  expect_no_error(visualize_result_box(result$epidish$rpc))
  # wrong input
  expect_no_error(visualize_result_bar(0))
  # correct result
  expect_no_error(visualize_result_bar(result$epidish$cbs))
  # wrong input
  expect_no_error(compare_results(0,result$tca))
  # correct result
  expect_no_error(compare_results(result$flowsortedbloodepic, result$methylcc))
})
