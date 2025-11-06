# Converts azimuth angles into Cardinal directions

Converts azimuth angles into Cardinal directions

## Usage

``` r
azimuth_to_cardinal(x, n_directions = 8)
```

## Arguments

- x:

  angles in degree.

- n_directions:

  either 4, for 4-points (N, E, S, W), 8 for 8-point (N, NE, E, …) or 6
  for 16-point (N, NNE, NE, …) cardinal version.

## Value

character vector

## See also

Other parse-orientations:
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md),
[`quadrant2dd()`](https://tobiste.github.io/structr/reference/quadrant2dd.md),
[`sense_from_strikeslip()`](https://tobiste.github.io/structr/reference/sense_from_strikeslip.md),
[`split_trailing_letters()`](https://tobiste.github.io/structr/reference/split_trailing_letters.md),
[`strikeslip_kinematics()`](https://tobiste.github.io/structr/reference/strikeslip_kinematics.md)

## Examples

``` r
azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 4) # 8-point compass
#> [1] "N" "N" "E" "E" "S" "W" "N"
azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 8) # 8-point compass
#> [1] "N"  "NE" "NE" "E"  "S"  "W"  "N" 
azimuth_to_cardinal(c(0, 23, 45, 100, 190, 270, 350), 16) # 16-point compass
#> [1] "N"   "NNE" "NE"  "E"   "S"   "W"   "N"  
```
