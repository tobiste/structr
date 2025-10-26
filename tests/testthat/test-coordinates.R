test_that("Line to cartesian works", {
  expect_equal(Vec3(Line(90, 0)), Vec3(0, 1, 0))
})

test_that("Cartesian to Line works", {
  expect_equal(Line(Vec3(1, 0, 0)), Line(0, 0))
  expect_equal(Line(Vec3(3, 2, 1)), Line(33.69007, 15.50136), tolerance = 1e-5)
})

