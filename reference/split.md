# Parse measurement and direction strings

Parse measurement and direction strings

## Usage

``` r
split_trailing_letters(x)

split_strike(x)
```

## Arguments

- x:

  character or number

## Value

list

## See also

Other parse-orientations:
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md),
[`azimuth_to_cardinal()`](https://tobiste.github.io/structr/reference/azimuth_to_cardinal.md),
[`quadrant2dd()`](https://tobiste.github.io/structr/reference/quadrant2dd.md),
[`sense_from_strikeslip()`](https://tobiste.github.io/structr/reference/sense_from_strikeslip.md),
[`strikeslip_kinematics()`](https://tobiste.github.io/structr/reference/strikeslip_kinematics.md)

## Examples

``` r
test <- c("45NW", "4SE")
split_strike(test)
#> $measurement
#> 45NW  4SE 
#>   45    4 
#> 
#> $direction
#> 45NW  4SE 
#> "NW" "SE" 
#> 
```
