test_that("Drillcore orientation works", {
  expect_equal(c(drillcore_orientation(225, -45, 60, 320)), c(25.003920, 70.029590), tolerance = 1e-6)
})
