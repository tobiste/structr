# Converts strike and dip quadrant notation into into dip direction

Converts strike and dip quadrant notation into into dip direction

## Usage

``` r
quadrant2dd(strike, dip_quadrant, n_directions = c(4L, 8L, 16L))
```

## Arguments

- strike:

  numeric. Strike in degree (left- or right-hand rule)

- dip_quadrant:

  character. Quadrant of dip direction

- n_directions:

  integer.

## Value

Dip direction in degrees

## See also

Other parse-orientations:
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md),
[`azimuth_to_cardinal()`](https://tobiste.github.io/structr/reference/azimuth_to_cardinal.md),
[`sense_from_strikeslip()`](https://tobiste.github.io/structr/reference/sense_from_strikeslip.md),
[`split()`](https://tobiste.github.io/structr/reference/split.md),
[`strikeslip_kinematics()`](https://tobiste.github.io/structr/reference/strikeslip_kinematics.md)

## Examples

``` r
s <- c(270, 315, 0, 45, 90, 135, 180, 225, 270) # strike in left-hand-rule
q <- c("N", "E", "E", "S", "S", "W", "W", "N", "N") # dip quadrant
quadrant2dd(s, q)
#> [1]   0  45  90 135 180 225 270 315   0
```
