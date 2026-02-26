# Projected Mean

Eigenvector with largest eigenvalue of the orientation tensor

## Usage

``` r
projected_mean(x, ...)
```

## Arguments

- x:

  either an object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`,
  `"Pair"`, or `"Fault"` where the rows are the observations and the
  columns are the coordinates, or an `"ortensor"` object.

- ...:

  additional arguments passed to
  [`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md)
  (ignored if `x` is `"ortensor"` object).

## Value

Vector in coordinate system of `x`

## References

Bachmann, F., Hielscher, R., Jupp, P. E., Pantleon, W., Schaeben, H., &
Wegert, E. (2010). Inferential statistics of electron backscatter
diffraction data from within individual crystalline grains. Journal of
Applied Crystallography, 43(6), 1338–1355.
https://doi.org/10.1107/S002188981003027X

## See also

[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)
for eigenvalues of orientation tensor,
[`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md) for
arithmetic mean,
[`geodesic_mean()`](https://tobiste.github.io/structr/reference/geodesic-mean.md)
for geodesic mean.

## Examples

``` r
example_lines_df$quality
#>  [1]  3  3 NA NA NA NA  4  4  4  3  3  3  3  3  3  3  3  3  3  3  3  3 NA  5  3
#> [26] NA  5  5  5  5  5  5  5  5  5  3  3  3  5  5  5  5  5  5  5  5  5  5  5  5
#> [51]  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
#> [76]  5  5  5  5  5  3  3  3  5
projected_mean(example_lines, w = runif(nrow(example_lines)))
#> Line object (n = 1):
#>  azimuth   plunge 
#> 67.80150 14.85824 
```
