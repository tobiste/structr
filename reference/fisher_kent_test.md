# Score test of the Fisher against the Kent distribution

Score test of the Fisher against the Kent distribution

## Usage

``` r
fisher_kent_test(x, level = 0.95)
```

## Source

modified after the python package `apsg` (v.2.0.6) by O. Lexa

## Arguments

- x:

  object of class `"Spherical"`

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

## Value

object of class `"Test"`

## Details

\\T = n (\kappa/2)^2 I_1/2(\kappa) / I_5/2(\kappa) r_2^2\\ is
asymptotically chi-square with 2 degrees of freedom (Kent, 1982, eq.
5.5), with the maximum-likelihood \\\kappa\\ of the Fisher distribution
(Kent's Theorem 5.1) and \\r_2\\ the difference of the two eigenvalues
of the second moments of the data in the tangent plane of the mean
direction.

Large \\T\\ indicates elliptical (Kent) rather than circular (Fisher)
contours.

## Note

FLE (sec. 5.4.6) computes the statistic with the concentration of the
fitted Kent distribution, which gives a larger value; this follows
Kent's theorem.

## References

Kent, J.T. (1982). The Fisher-Bingham distribution on the sphere.
Journal of the Royal Statistical Society B 44, 71-80.

## See also

Other tests:
[`bingham_test()`](https://tobiste.github.io/structr/reference/bingham_test.md),
[`fisher_ftest()`](https://tobiste.github.io/structr/reference/fisher_ftest.md),
[`fisher_goodness_of_fit()`](https://tobiste.github.io/structr/reference/fisher_goodness_of_fit.md),
[`rayleigh_test3()`](https://tobiste.github.io/structr/reference/rayleigh_test3.md),
[`watson_symmetry_test()`](https://tobiste.github.io/structr/reference/watson_symmetry_test.md),
[`wellner()`](https://tobiste.github.io/structr/reference/wellner.md)

## Examples

``` r
#' set.seed(20250411)
x <- rvmf(100)
fisher_kent_test(x)
#> 
#>  Fisher vs Kent
#> 
#> statistic = 0.12966, df = 2, p-value = 0.9372
#> critical value (95% confidence level) = 5.9915
#> null hypothesis rejected at the 95% confidence level: no
#> further statistics:
#>      kappa         r2 
#> 5.66200190 0.00954967 
#> source: Kent (1982), eq. 5.5 
#> 
```
