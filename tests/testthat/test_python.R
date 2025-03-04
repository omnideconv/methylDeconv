library(methylDeconv)

methyldeconv::init_python()

test_that("Python environment creation works", {
  
  expect_identical(info = "Python available", 
                   object = reticulate::py_available(), 
                   expected = TRUE)
})


test_that("Python environment exists", {

  expect_identical(info = "r-methyldeconv exists", 
                   object = reticulate::condaenv_exists('r-methyldeconv'), 
                   expected = TRUE)
})
