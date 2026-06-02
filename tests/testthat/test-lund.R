test_that("Lund works", {
  S1 <- Line(250.89, 70.07)
  S3 <- Line(103.01, 17.07)
  S2 <- crossprod(S3, S1)
  sh <- SH(S1, S2, S3, R = 1)
  
  expect_type(sh, "double")
  expect_equal(sh, 70.89)
})
