# Mean strain ellipse

determines the shape and orientation of the strain ellipse by using
deformed elliptical objects as strain markers. The algorithm is based on
the mean shape matrix and its eigenvalues (Shimamoto and Ikeda, 1976).

## Usage

``` r
mean_strain_ellipse(r, phi, boot = TRUE, resamples = 1000, boot.values = FALSE)
```

## Arguments

- r:

  numeric. Aspect ratio of deformed object (long axis / short axis)

- phi:

  numeric. Orientation of long axis of deformed object (in degrees)

- boot:

  logical. Whether a 95% confidence interval from on bootstrapping
  should be calculated. `TRUE` by default.

- resamples:

  integer. Number of bootstrap resamples (`1000` by default). Ignored
  when `boot = FALSE`.

- boot.values:

  logical. Whether the bootstrapped R and phi values should be added to
  the output. `FALSE` by default.

## Value

list. `R` gives the mean aspect ratio of the strain ellipse, and `phi`
gives the orientation of its long axis. If `boot=TRUE`, then the
bootstrapped 95% confidence interval for the mean aspect ratio (`R_CI`)
and for its orientation (`phi_CI`) are added. If `boot.values=TRUE`, a
matrix containing the bootstrapped R and phi values are added.

## References

Shimamoto, T., Ikeda, Y., 1976. A simple algebraic method for strain
estimation from ellipsoidal objects. Tectonophysics 36, 315–337.
[doi:10.1016/0040-1951(76)90107-4](https://doi.org/10.1016/0040-1951%2876%2990107-4)

## See also

[`mean_strain_ellipse_ramsay()`](https://tobiste.github.io/structr/reference/mean_strain_ellipse_ramsay.md)

## Examples

``` r
set.seed(20250411)
data(ramsay)
mean_strain_ellipse(ramsay[, "R"], ramsay[, "phi"])
#> $R
#> [1] 1.628138
#> 
#> $phi
#> [1] 25.73632
#> 
#> $R_CI
#> [1] 1.59275 1.66359
#> 
#> $phi_CI
#> [1] 24.71905 26.75180
#> 
```
