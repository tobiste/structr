# excel
test_that('import excel', {
  expect_no_failure(read_strabo_xls('strabo_example.xlsx'))
})

# mobile
test_that('import csv', {
  expect_no_failure(read_strabo_mobile('strabo_example.txt'))
})


# json
test_that('import json', {
  expect_no_failure(read_strabo_JSON('strabo_example.json'))
})
