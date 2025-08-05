f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))

test_that("multiplication works", {
  expect_equal(Fault_rake(f), c(84.72020, -10.28562, 30.11825), tolerance = 1e-5)
})
