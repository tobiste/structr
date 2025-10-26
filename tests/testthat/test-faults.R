f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))

test_that("multiplication works", {
  expect_equal(Fault_rake(f), c(84.72020, -10.28562, 30.11825), tolerance = 1e-5)
})




# test from Shorrot and Lise, 1998, Tab. 2
strike_con7 <- c(270, 315, 0, 45, 90, 135, 180, 225, 270) # strike in left-hand-rule
strike_con6 <- c(90,135,180,225,270,315,0,45,90)
quadrant_con4 <- c("N", 'E', 'E', 'S', 'S', 'W', 'W', 'N', "N") # dip quadrant
dip_direction_con2 <- seq(0, 360, 45) %% 360

test_that("Quadrant convenction", {
  expect_equal(quadrant2dd(strike_con7, quadrant_con4), dip_direction_con2)
  expect_equal(quadrant2dd(strike_con6, quadrant_con4), dip_direction_con2)
  
})


strike <- c(90, 135, 180, 225, 270, 270, 315, 0, 45)
dip_dir <- c(180, 225, 270, 315, 360, 0, 45, 90, 135) %% 360
test_that("Right-hand rule 2", {
  expect_equal(rhr2dd(strike), dip_dir)
})

dip_quadrant <- c("S", 'S', 'W', 'W', 'N', 'N', 'N', 'E', 'E')
# test_that("dip quadrant", {
#   expect_equal(azimuth_to_cardinal(dip_dir, 4), dip_quadrant)
# })


dip <- c(5, 10, 15, 30, 40, 55, 65, 75, 89.999)
rake <- seq(0, 360, 45) %% 360

f <- Fault_from_rake(Plane(dip_dir, dip), rake)

azimuth <- c(90, 179.6, 270, 4.1, 90, 299.8, 45, 165.5, 45)
test_that("Test Fault from rake: azimuth", {
  expect_equal(round(f[, 3], 1), azimuth)
})

plunge <- c(0, 7.1, 15, 20.7, 0, 144.6, 115, 115, 136.9, 180)
# test_that("Test Fault from rake: plunge", {
#   expect_equal(round(Line(f)[, 2] %% 360, 1), plunge)
# })






# test from Shorrot and Lise, 1998, Tab. 2
# dip_con1 <- c(5, 10, 15, 30, 40, 55, 65, 75, 90)
# dip_dir_con2 <- c(180, 225, 270, 315, 360, 0, 45, 90, 135)
# 
# 
# plunge_11 <- c(0, 7, 15, 21, 0, 35, 65, 43, 0)
# azimuth_12 <- c(90, 179, 270, 4, 90, 299, 45, 165, 225)
# 
# pitch_con17 <- c(0, 45, 90, 135, 180, 45, 90, 135, 180)
# pitch_con15 <- c(0, 45, 90, 45, 0, 45, 90, 45, 0)
# plunge_quadrant_14 <- c("E", "S", "W", "N", "E", "W", "E", "S", "W")
# pitch_quadrant_con16 <- c("E", "S", "S", "E", "E", 'W', "N", "S", "W")
# rake_con_18 <- dip_direction_con2
# exp_sense <- ifelse(rake_con_18 <= 180, 1, -1)

# test_that("Quadrant convenction rake 1", {
#   expect_equal(
#     round(Fault_from_rake_quadrant(Plane(dip_dir_con2, dip_con1), pitch_con17, plunge_quadrant_14, "plunge")|> Fault_rake()), 
#     rake_con_18)
# })
# 
# test_that("Quadrant convenction rake 2", {
#   expect_equal(
#     round(Fault_from_rake_quadrant(Plane(dip_dir_con2, dip_con1), pitch_con15, pitch_quadrant_con16, "rake") |> Fault_rake()), 
#     rake_con_18)
# })
# 
# test_that("Fault rake", {
#   expect_equal(
#     round(Fault_rake(Fault(dip_dir_con2, dip_con1, azimuth_12, plunge_11, sense = exp_sense)) %% 360),
#     rake_con_18,
#     tolerance = 2/360
#   )
# })
