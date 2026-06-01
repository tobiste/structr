# excel
test_that("import excel", {
  expect_s3_class(read_strabo_xls(test_path("testdata", "strabo_example.xlsx")), 'strabo')
})

# mobile
test_that("import csv", {
  expect_s3_class(read_strabo_mobile(test_path("testdata", "strabo_example.txt")), 'strabo')
})

# json
test_that("import json", {
  expect_s3_class(read_strabo_JSON(test_path("testdata", "strabo_example.json")), 'strabo')
})
