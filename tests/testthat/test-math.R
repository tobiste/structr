test_that("Spherical mean", {
  expect_equal(c(v_mean(Line(c(0, 0), c(0, 90)))), c(0, 45))
})

test_that("Spherical variance", {
  expect_equal(v_var(Line(c(0, 0), c(0, 90))), 0.29289322)
})

test_that("Spherical standard deviation", {
  expect_equal(v_delta(Line(c(0, 0), c(0, 90))), 45)
})
