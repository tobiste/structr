# Ellipsoid tensor from principal stretches

Return diagonal tensor defined by magnitudes of principal stretches

## Usage

``` r
ellipsoid_from_stretch(x = 1, y = 1, z = 1)
```

## Arguments

- x, y, z:

  numeric. Magnitudes of principal stretches

## Value

object of class `"ellipsoid"`

## Examples

``` r
el <- ellipsoid_from_stretch(4, 3, 1)
principal_stretch(el)
#> S1 S2 S3 
#>  4  3  1 
```
