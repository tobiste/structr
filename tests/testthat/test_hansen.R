# ds6 <- read.delim("C:/Users/tstephan/Downloads/hanson_ds6.txt", sep = '\t', header = F) |> 
#   mutate(V1 = rhr2dd(V1)) |> 
#   as.Pair()
# 
# slip_inversion_hansen(ds6, 1)
test_that("import json", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.json')
  data_path <- test_path("testdata", "hanson_ds6.txt")
  
  import <- as.matrix(read.delim(data_path, sep = '\t', header = FALSE))
  import[, 1] <- rhr2dd(import[, 1])
  
  res <- slip_inversion_hansen(as.Pair(import), TRUE)
  
  # expect_s3_class(res, 'list')
  expect_equal(round(res$stress_shape$phi, 2), 0.18)
  expect_equal(round(res$vorticity_mag, 2), -0.8)
  expect_equal(c(round(res$principal_axes[1, ])), c(24, 87))
  expect_equal(c(round(res$principal_axes[2, ])), c(244, 2))
  expect_equal(c(round(res$principal_axes[3, ])), c(154, 2))
  expect_equal(c(round(res$vorticity_axis)), c(240, 3))
})
