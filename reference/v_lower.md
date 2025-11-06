# Antipodal map sending upper-hemisphere rays to the lower hemisphere.

Negates any vector with positive `z`-component. Used in lower-hemisphere
plots.

## Usage

``` r
v_lower(v)
```

## Arguments

- v:

  object of class `'Vec3'`, where the rows are the observations and the
  columns are the Cartesian coordinates.

## Value

object of class `'Vec3'`

## Examples

``` r
v <- rbind(Vec3(1, 0, -1), Vec3(1, 0, 1))
v_lower(v)
#> Vector (Vec3) object (n = 2):
#>       x y  z
#> [1,]  1 0 -1
#> [2,] -1 0 -1
```
