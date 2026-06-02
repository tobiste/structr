# excel
test_that("import excel", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.xlsx')
  data_path <- test_path("testdata", "strabo_example.xlsx")
  
  import <- read_strabo_xlsx(data_path)
  expect_s3_class(import, 'strabo')
  expect_equal(nrow(import$data), 3)
  expect_equal(nrow(import$linear), 3)
  expect_equal(nrow(import$planar), 3)
})

test_that("subset strabo excel", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.xlsx')
  data_path <- test_path("testdata", "strabo_example.xlsx")
  
  import <- read_strabo_xlsx(data_path)
  import_sub <- subset(import, !is.na(Orientation.Unix.Timestamp))
  
  expect_s3_class(import_sub, 'strabo')
  expect_equal(nrow(import_sub$data), 2)
  expect_equal(nrow(import_sub$linear), 2)
  expect_equal(nrow(import_sub$planar), 2)
})


# mobile
test_that("import txt", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.txt')
  data_path <- test_path("testdata", "strabo_example.txt")
  
  import <- read_strabo_mobile(data_path)
  expect_s3_class(import, 'strabo')
  expect_equal(nrow(import$data), 4)
  expect_equal(nrow(import$linear), 4)
  expect_equal(nrow(import$planar), 4)
})

test_that("subset strabo txt", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.txt')
  data_path <- test_path("testdata", "strabo_example.txt")
  
  import <- read_strabo_mobile(data_path)
  import_sub <- subset(import, Notes == 'a note')
  
  expect_s3_class(import_sub, 'strabo')
  expect_equal(nrow(import_sub$data), 1)
  expect_equal(nrow(import_sub$linear), 1)
  expect_equal(nrow(import_sub$planar), 1)
})


# json
test_that("import json", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.json')
  data_path <- test_path("testdata", "strabo_example.json")
  
  import <- read_strabo_JSON(data_path)

  expect_s3_class(import, 'strabo')
  expect_equal(nrow(import$data), 23)
  expect_equal(nrow(import$linear), 23)
  expect_equal(nrow(import$planar), 23)                
})


test_that("subset strabo json", {
  # data_path <- file.path('tests', 'testthat', 'testdata', 'strabo_example.json')
  data_path <- test_path("testdata", "strabo_example.json")
  import <- read_strabo_JSON(data_path)
  import_sub <- subset(import, feature_type == 'bedding')
  
  expect_s3_class(import_sub, 'strabo')
  expect_equal(nrow(import_sub$data), 11)
  expect_equal(nrow(import_sub$linear), 11)
  expect_equal(nrow(import_sub$planar), 11)                
})