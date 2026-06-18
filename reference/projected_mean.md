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
projected_mean(example_lines, w = runif(nrow(example_lines)))
#> Line object (n = 1):
#>  azimuth   plunge 
#> 68.41056 13.42329 
```
