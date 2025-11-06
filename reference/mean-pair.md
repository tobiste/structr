# Mean orientation of a set of pairs or faults

The Frechet (geodesic \\L^2\\) mean and variance of a pair of foliations
and lineations

## Usage

``` r
geodesic_mean_pair(x, group = NULL)

geodesic_var_pair(x, group = NULL)
```

## Source

oriMeanVariance from geologyGeometry (J.R. Davis)

## Arguments

- x:

  object of class `"Pair"` or `"Fault"`

- group:

  Symmetry group of `x`. See
  [`symmetry_group()`](https://tobiste.github.io/structr/reference/symmetry_group.md)
  for details If `NULL`, the group will be automatically picked based on
  the class of `x`.

## Value

object of class `"Pair"` or `"Fault"`, respectively

## References

Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis for
three-dimensional orientational data. Journal of Structural Geology, 96,
65–89. https://doi.org/10.1016/j.jsg.2017.01.002

## See also

[`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md) for
arithmetic mean, and
[`geodesic_mean_line()`](https://tobiste.github.io/structr/reference/geodesic-line.md)
for geodesic mean of lines.

## Examples

``` r
my_fault <- Fault(
  c("a" = 120, "b" = 120, "c" = 100),
  c(60, 60, 50),
  c(110, 25, 30),
  c(58, 9, 23),
  c(1, -1, 1)
)
geodesic_mean_pair(my_fault)
#> Fault object (n = 1):
#> dip_direction           dip       azimuth        plunge         sense 
#>     109.11409      59.36207     131.87199      57.28751       1.00000 
geodesic_var_pair(my_fault)
#> [1] 0.721249
```
