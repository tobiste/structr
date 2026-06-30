test_that("import json", {
  res <- slip_inversion_hansen(osmundsen2010, TRUE)
  
  # expect_s3_class(res, 'list')
  expect_equal(round(res$stress_shape$phi, 2), 0.18)
  expect_equal(round(res$vorticity_mag, 2), -0.8)
  expect_equal(c(round(res$principal_axes[1, ])), c(24, 87))
  expect_equal(c(round(res$principal_axes[2, ])), c(244, 2))
  expect_equal(c(round(res$principal_axes[3, ])), c(154, 2))
  expect_equal(c(round(Line(res$vorticity_axis))), c(240, 3))
})
