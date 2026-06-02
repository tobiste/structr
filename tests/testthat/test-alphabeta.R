test_that("Drillcore orientation works", {
  dat <- drillcore_transformation(225, 45, 60, 320)
  expect_s3_class(dat, "Plane")
  expect_equal(c(dat), c(25.003920, 70.029590), tolerance = 1e-6)
})
