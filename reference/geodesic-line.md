# The Fréchet (geodesic \\L^2\\) mean of a set of lines or rays

An iterative algorithm for computing the Fréchet mean — the line or ray
that minimizes the Fréchet variance. The iterations continue until error
squared of epsilon is achieved or `steps` iterations have been used. Try
multiple seeds, to improve your chances of finding the global optimum.

## Usage

``` r
geodesic_mean_line(x, ...)

geodesic_mean_ray(x, ...)

geodesic_var_line(x, ...)

geodesic_var_ray(x, ...)

geodesic_meanvariance_line(x, seeds = 5L, steps = 100L)

geodesic_meanvariance_ray(x, seeds = 5L, steps = 100L)
```

## Source

lineMeanVariance from geologyGeometry (J.R. Davis)

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- ...:

  parameters passed to `geodesic_meanvariance_line()` or
  `geodesic_meanvariance_ray()`

- seeds:

  positive integer. How many `x` to try as seeds

- steps:

  positive integer. Bound on how many iterations to use.

## Value

`geodesic_meanvariance_*` (\* either line or ray) returns a `list`
consisting of `$variance` (numeric), `$mean` (a line), `$error` (an
integer) and `$min.eigenvalue` (numeric). `geodesic_mean_*` and
`geodesic_var_*` are convenience wrapper and only return the mean and
the variance, respectively.

## Details

Error should be `0` and `min.eigenvalue` should be positive. Otherwise
there was some problem in the optimization. If error is non-zero, then
try increasing `steps`.

## References

Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis for
three-dimensional orientational data. Journal of Structural Geology, 96,
65–89. https://doi.org/10.1016/j.jsg.2017.01.002

## See also

[`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md) for
arithmetic mean and
[`geodesic_mean_pair()`](https://tobiste.github.io/structr/reference/mean-pair.md)
for geodesic mean of pairs.

## Examples

``` r
geodesic_meanvariance_line(example_lines)
#> $variance
#> [1] 0.06118261
#> 
#> $mean
#> Line object (n = 1):
#>  azimuth   plunge 
#> 69.64251 14.87961 
#> 
#> $error
#> [1] 0
#> 
#> $min.eigenvalue
#> [1] 0
#> 
geodesic_mean_line(example_lines)
#> Line object (n = 1):
#>  azimuth   plunge 
#> 69.64075 14.87885 
geodesic_var_line(example_lines)
#> [1] 0.06118261
```
