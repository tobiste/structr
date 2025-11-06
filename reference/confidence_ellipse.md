# Confidence ellipse

Bootstrapped projected mean with percentile confidence region and
hypothesis tests

## Usage

``` r
confidence_ellipse(x, n = 10000L, alpha = 0.05, res = 512L, isotropic = FALSE)
```

## Source

geologyGeometry (J.R. Davis)

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- n:

  integer (10000 by default).

- alpha:

  numeric (0.05 by default).

- res:

  integer. The resolution with which to sample the boundary curve of the
  (1 - `alpha`) \* 100% percentile region.

- isotropic:

  logical. If `TRUE`, forces the inverse covariance to the identity
  matrix and hence the region to be circular.

## Value

list.

- `center`:

  Projected mean of `x`

- `cov`:

  Inverse covariance matrix in the tangent space, which is the identity
  if `isotropic` is `TRUE`

- `rotation`:

  Rotation matrix used

- `quantiles`:

  Quantiles of Mahalanobis norm

- `pvalue.ray`:

  For rays: p-value for each ray in `x`, i.e. the fraction of `x` that
  are farther from `center` than the given ray.

- `pvalue.line`:

  For lines: p-value for each line in `x`, i.e. the fraction of `x` that
  are farther from `center` than the given line

- `pvalue.line.FUN`,`pvalue.ray.FUN`:

  The function to calculate the p-value for a given vector

- `angles`:

  Angles of the semi-axis of the confidence ellipse (in radians if `x`
  is an `"Vec3"` object, in degrees if otherwise.)

- `ellipse`:

  Confidence ellipse given as `"Vec3"` object with `res` vectors

## Note

The inference is based on percentiles of Mahalanobis distance in the
tangent space at the mean of the bootstrapped means. The user should
check that the bootstrapped means form a tight ellipsoidal cluster,
before taking such a region seriously.

## References

Davis, J. R., & Titus, S. J. (2017). Modern methods of analysis for
three-dimensional orientational data. Journal of Structural Geology, 96,
65–89.
[doi:10.1016/j.jsg.2017.01.002](https://doi.org/10.1016/j.jsg.2017.01.002)

## Examples

``` r
set.seed(20250411)
ce <- confidence_ellipse(example_lines, n = 1000, res = 10)
# print(ce)

# Check how many vectors lie outside quantiles:
stats::quantile(ce$pvalue, probs = c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00))
#>      0%      5%     25%     50%     75%    100% 
#> 0.00000 0.04995 0.24975 0.49950 0.74925 0.99900 

# Hypothesis testing (reject if p-value < alpha):
ce$pvalue.FUN((Line(90, 0)))
#> [1] 0
```
