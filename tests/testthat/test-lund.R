S1 <- Line(250.89, 70.07)
S3 <- Line(103.01, 17.07)
S2 <- vcross(S3, S1)
SH(S1, S2, S3, R = 1)

test_that("Lund works", {
  expect_equal(SH(S1, S2, S3, R = 1), 70.89)
})
