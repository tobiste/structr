# Dip-slip Kinematics from Strike-Slip Faults

Dip-slip sense (1 for normal, -1 for reverse) when only strike-slip
kinematics are known

## Usage

``` r
sense_from_strikeslip(x, left)
```

## Arguments

- x:

  `"Pair"` object(s) of the plane(s) and line(s) with unknown dip-slip
  offset

- left:

  logical. `TRUE` if `x` is a left-lateral (sinistral) fault, and
  `FALSE` if `x` is a right-lateral (dextral) fault. Must have same
  length as rows in `x`

## Value

numeric. 1 if normal, -1 if reverse offset

## See also

[`strikeslip_kinematics()`](https://tobiste.github.io/structr/reference/strikeslip_kinematics.md)

Other parse-orientations:
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md),
[`azimuth_to_cardinal()`](https://tobiste.github.io/structr/reference/azimuth_to_cardinal.md),
[`quadrant2dd()`](https://tobiste.github.io/structr/reference/quadrant2dd.md),
[`split_trailing_letters()`](https://tobiste.github.io/structr/reference/split_trailing_letters.md),
[`strikeslip_kinematics()`](https://tobiste.github.io/structr/reference/strikeslip_kinematics.md)

## Examples

``` r
# Sinistral fault
sense_from_strikeslip(Pair(120, 89, 30, 5), left = TRUE) # 1: normal offset
#> [1] 1

# Dextral fault
sense_from_strikeslip(Pair(120, 89, 30, 5), left = FALSE) # -1: reverse offset
#> [1] -1
```
