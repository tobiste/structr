test_that("Line to cartesian works", {
  expect_equal(c(line2vec(Line(90, 0))), c(0, 1, 0))
})

test_that("Cartesian to Line works", {
  expect_equal(c(vec2line(c(1, 0, 0))), c(0, 0))
  expect_equal(c(vec2line(c(3, 2, 1))), c(33.69007, 15.50136), tolerance = 1e-5)
})
