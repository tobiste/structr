# The Fréchet (geodesic \\L^2\\) mean

An iterative algorithm for computing the Fréchet mean, i.e. the vector
that minimizes the Fréchet variance.

## Usage

``` r
geodesic_mean(x, ...)

# S3 method for class 'Vec3'
geodesic_mean(x, ...)

# S3 method for class 'Line'
geodesic_mean(x, ...)

# S3 method for class 'Ray'
geodesic_mean(x, ...)

# S3 method for class 'Plane'
geodesic_mean(x, ...)

# S3 method for class 'Pair'
geodesic_mean(x, ...)
```

## Source

geologyGeometry (J.R. Davis)

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- ...:

  parameters passed to
  [`geodesic_meanvariance_ray()`](https://tobiste.github.io/structr/reference/geodesic-line.md)
  (if `x` is a Ray),
  [`geodesic_meanvariance_line()`](https://tobiste.github.io/structr/reference/geodesic-line.md)
  (if `x` is a Vec3, Line or Plane) or
  [`geodesic_mean_pair()`](https://tobiste.github.io/structr/reference/mean-pair.md)
  (if `x` is a Pair or a Fault).

## Value

the Fréchet mean vector as an object of class `x`.

## References

Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis for
three-dimensional orientational data. Journal of Structural Geology, 96,
65–89.
[doi:10.1016/j.jsg.2017.01.002](https://doi.org/10.1016/j.jsg.2017.01.002)

## See also

[`geodesic_var()`](https://tobiste.github.io/structr/reference/geodesic-var.md)
for Fréchet variance,
[`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md) for
the arithmetic mean,
[`projected_mean()`](https://tobiste.github.io/structr/reference/projected_mean.md)
for projected mean

## Examples

``` r
set.seed(20250411)
geodesic_mean(example_planes)
#> Plane object (n = 1):
#> dip_direction           dip 
#>     345.92363      75.35166 
```
