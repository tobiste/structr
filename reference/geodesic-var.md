# The Fréchet (geodesic \\L^2\\) variance

Dispersion measured using the Fréchet variance, i.e the sum of the
squared geodesic distances between all vectors and a specified vector.

## Usage

``` r
geodesic_var(x, ...)

# S3 method for class 'Vec3'
geodesic_var(x, y = NULL, ...)

# S3 method for class 'Line'
geodesic_var(x, y = NULL, ...)

# S3 method for class 'Plane'
geodesic_var(x, y = NULL, ...)

# S3 method for class 'Ray'
geodesic_var(x, y = NULL, ...)

# S3 method for class 'Pair'
geodesic_var(x, y = NULL, group = NULL, ...)
```

## Source

`geologyGeometry` (J.R. Davis): http://www.joshuadavis.us/software/

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

- y:

  Only for variance. object of class `"Vec3"`, `"Line"`, `"Ray"`,
  `"Plane"`, `"Pair"`, or `"Fault"` about which the Fréchet variance
  should be calculated for. If `NULL` (the default), Fréchet variance
  about the Fréchet mean.

- group:

  Symmetry group of `x`. See
  [`symmetry_group()`](https://tobiste.github.io/structr/reference/symmetry_group.md)
  for details If `NULL`, the group will be automatically picked based on
  the class of `x`.

## Value

the Fréchet variance as a numeric number. Because distances in SO(3)
never exceed \\\pi\\, the maximum possible variance is \\\frac{\pi^2}{2}
\approx 4.93\\.

## Details

The variance of a dataset \\{x_1, \ldots, x_n}\\ about a vector \\y\\ is
defined as \$\$ \Psi(x) = \frac{1}{2n} \sum\_{i=1}^n d_G(y, x_i)^2\$\$
where \\d_G(x, y)\\ is the geodesic distance between vectors \\x\\ and
\\y\\ (see
[`angle()`](https://tobiste.github.io/structr/reference/vecmath.md)).

## References

Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis for
three-dimensional orientational data. Journal of Structural Geology, 96,
65–89.
[doi:10.1016/j.jsg.2017.01.002](https://doi.org/10.1016/j.jsg.2017.01.002)

## See also

[`geodesic_mean()`](https://tobiste.github.io/structr/reference/geodesic-mean.md)
for the Fréchet mean,
[`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md) for
the arithmetic mean,
[`projected_mean()`](https://tobiste.github.io/structr/reference/projected_mean.md)
for projected mean

## Examples

``` r
set.seed(20250411)
geodesic_var(example_planes, example_planes[1, ])
#> [1] 0.559469
geodesic_var(example_planes)
#> [1] 0.2656372
```
