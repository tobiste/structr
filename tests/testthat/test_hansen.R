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
  
  res <- slip_inversion_hansen(as.Pair(import), FALSE)
  
  # expect_s3_class(res, 'list')
  expect_equal(round(res$nine$phi, 2), 0.18)
  expect_equal(round(res$six$phi, 2), 1-0.02)
  # expect_equal(round(res$nine$vorticity_mag, 2), -0.8)
})
