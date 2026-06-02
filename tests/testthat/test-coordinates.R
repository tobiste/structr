test_that("objects class return", {
  expect_s3_class(Vec3(1, 0, 0), "Vec3")
  expect_s3_class(Line(90, 0), "Line")
  expect_s3_class(Plane(90, 0), "Plane")
  expect_s3_class(Pair(90, 0, 90, 0), "Pair")
  expect_s3_class(Fault(90, 0, 90, 0, 1), "Fault")
  expect_s3_class(Fault_from_rake(Plane(90, 0), 12, 1), "Fault")
})


test_that("Line to cartesian works", {
  expect_equal(Vec3(Line(90, 0)), Vec3(0, 1, 0))
})

test_that("Cartesian to Line works", {
  expect_equal(Line(Vec3(1, 0, 0)), Line(0, 0))
  expect_equal(Line(Vec3(3, 2, 1)), Line(33.69007, 15.50136), tolerance = 1e-5)
})
