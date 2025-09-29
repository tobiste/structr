## code to prepare `DATASET` dataset goes here
# simongomez = read.delim("clipboard")
# simongomez <- simongomez |> mutate(
#   dipdir = Strike + 90,
#   ve = NA, ve = ifelse(sense %in% c("D", "T"), 1, ve), ve = ifelse(sense %in% c("S", "N"), -1, ve))
# fault = Plane(simongomez$dipdir, simongomez$Dip)
# theta = simongomez$Rake
# ve = simongomez$ve
# usethis::use_data(simongomez, overwrite = TRUE)

# readr::read_csv("Field Work/data_planes.txt") |>
#   dplyr::filter(cluster == "Huronian Lk.") |>
#   select(dipdir, dip, quality, feature_type) |>
#   write_csv("C:/Users/tstephan/Documents/GitHub/structr/inst/example_planes.csv")
#
# lines <- readr::read_csv("Field Work/data_lines.txt") |>
#   dplyr::filter(cluster == "Huronian Lk.") |>
#   select(trend, plunge, quality, feature_type) |>
#   write_csv("C:/Users/tstephan/Documents/GitHub/structr/inst/example_lines.csv")
example_planes <- readr::read_csv("inst/example_planes.csv")
usethis::use_data(example_planes, overwrite = TRUE)

example_lines <- readr::read_csv("inst/example_lines.csv")
usethis::use_data(example_lines, overwrite = TRUE)


ramsay <- read.csv("inst/RamsayHuber1983.csv", header = FALSE)
colnames(ramsay) <- c("R", "phi")
ramsay <- as.matrix(ramsay)
usethis::use_data(ramsay, overwrite = TRUE)



holst <- tibble::tribble(
  ~locality, ~e1e2, ~e2e3,
  1, 0.03, 3.43,
  2, 0.07, 3.3,
  3, 0.04, 3.22,
  4, 0.07, 3.33,
  5, 0.14, 3,
  6, 0.75, 2.2,
  7, 0.96, 1.61,
  8, 1.64, 1.25,
  9, 0.56, 2.08,
  10, 1.45, 1.39,
  11, 1.95, 1.1,
  12, 2.03, 0.92,
  13, 2.69, 0.83,
  14, 2.2, 1.38,
  15, 1.76, 1.13,
  16, 0.13, 2.64,
  17, 0.05, 2.89
) |>
  # mutate(locality = as.integer(locality)) |>
  mutate(
    R_XY = exp(e1e2),
    R_YZ = exp(e2e3),
    R_XZ = R_XY * R_YZ,
    K = e1e2 / e2e3, # Hossack 1968
    v = (1 - K) / (1 + K),
    es = 1 / sqrt(3) * sqrt(log(R_XY)^2 + log(R_YZ)^2 + log(1 / R_XZ)^2), # Nadai, 1963
    r = R_XY + R_YZ - 1 # Wattersion 1968
  ) |>
  select(R_XY, R_YZ) |>
  as.matrix()
# holst
usethis::use_data(holst, overwrite = TRUE)

# Compare with Holst and Fossen (1987), Fig. 5:
# plot(log(holst[, 'R_YZ']), log(holst[, 'R_XY']), asp = 1, xlim = c(0, 4), ylim = c(0, 4))
