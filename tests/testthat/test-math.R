test_that("Spherical mean", {
  l <- Line(c(0, 0), c(0, 90))
  lm <- sph_mean(l)
  
  expect_s3_class(lm, "Line")
  expect_equal(c(lm), c(0, 45))
})

test_that("Spherical variance", {
  l <- Line(c(0, 0), c(0, 90))
  lv <- sph_var(l)
  
  expect_type(lv, "double")
  expect_equal(lv, 0.29289322)
})

test_that("Spherical standard deviation", {
  l <- Line(c(0, 0), c(0, 90))
  d <- delta(l)
  
  expect_type(d, "double")
  expect_equal(d, 45)
})
