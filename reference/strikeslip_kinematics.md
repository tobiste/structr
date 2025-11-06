# Strike-slip Kinematics

Returns the strike-slip kinematics of a fault

## Usage

``` r
strikeslip_kinematics(x)
```

## Arguments

- x:

  `"Fault"` object(s)

## Value

character. `"left"` - left-lateral (sinistral) offset, `"right"` -
right-lateral (dextral) offset

## See also

[`sense_from_strikeslip()`](https://tobiste.github.io/structr/reference/sense_from_strikeslip.md)

Other parse-orientations:
[`Fault_from_rake_quadrant()`](https://tobiste.github.io/structr/reference/Fault_from_rake_quadrant.md),
[`azimuth_to_cardinal()`](https://tobiste.github.io/structr/reference/azimuth_to_cardinal.md),
[`quadrant2dd()`](https://tobiste.github.io/structr/reference/quadrant2dd.md),
[`sense_from_strikeslip()`](https://tobiste.github.io/structr/reference/sense_from_strikeslip.md),
[`split_trailing_letters()`](https://tobiste.github.io/structr/reference/split_trailing_letters.md)

## Examples

``` r
strikeslip_kinematics(Fault(120, 50, 104, 49, sense = -1))
#> [1] "right"
```
