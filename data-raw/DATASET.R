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
example_planes_df <- readr::read_csv("inst/example_planes.csv")
usethis::use_data(example_planes_df, overwrite = TRUE)

example_planes <- Plane(example_planes_df$dipdir, example_planes_df$dip)
usethis::use_data(example_planes, overwrite = TRUE)


example_lines_df <- readr::read_csv("inst/example_lines.csv")
usethis::use_data(example_lines_df, overwrite = TRUE)

example_lines <- Line(example_lines_df$trend, example_lines_df$plunge)
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
# 
# 
# 
# 
gr <- tibble::tribble(
  ~strike, ~dip, ~dipdir_quadrant, ~facing, ~type,
  142, 58, "W", "L", "Cleavage",
  155, 63, "W", "L", "Cleavage",
  166, 82, "W", "L", "Cleavage",
  177+180, 83, "E", "U", "Cleavage",
  9, 68, "E", "U", "Cleavage",
  11, 60, "E", "U", "Cleavage",
  17, 48, "E", "U", "Cleavage",
  150, 65, "W", "L", "Cleavage",
  156, 50, "W", "U", "Bedding",
  146, 45, "W", "U", "Bedding",
  138, 34, "W", "U", "Bedding",
  98, 33, "W", "U", "Bedding",
  62, 36, "E", "U", "Bedding",
  43, 48, "E", "U", "Bedding",
  34, 48, "E", "U", "Bedding",
  29, 61, "E", "U", "Bedding"
) |>
  mutate(dipdir = rhr2dd(strike), .after = dip) |> 
  # mutate(dipdir_quadrant2 = azimuth_to_cardinal(dipdir), .after = dipdir_quadrant) |> 
  mutate(sense = ifelse(facing=="L", 1, -1))

# gray_example <- Ray(gr$dipdir, gr$dip, sense = gr$sense) |> Line() |> as.Plane()
gray_example <- Plane(gr$dipdir, gr$dip, sense = gr$sense)
rownames(gray_example) <- gr$type
usethis::use_data(gray_example, overwrite = TRUE)



# Angelier 1990 fault slip data
angelier1990_df <- readxl::read_xlsx("inst/angelier_data.xlsx", sheet = 1) |>
  dplyr::group_by(site) |>
  dplyr::mutate(id = dplyr::row_number(), .before = "sense") |>
  dplyr::mutate( # dipdir = structr::rhr2dd(strike),
    sense = ifelse(sense == "N", 1, -1)
  )
angelier1990 <- angelier1990_df |> split(angelier1990_df$site) |> 
  lapply(function(x){
    Fault_from_rake_quadrant(
      Plane(quadrant2dd(x$strike, x$dip_card), x$dip), 
      rake = x$pitch, 
      quadrant = x$pitch_card, 
      type = "rake", 
      sense = x$sense)
})
usethis::use_data(angelier1990, overwrite = TRUE)


