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
