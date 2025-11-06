# Mean strain ellipse after Ramsay (1967)

Mean strain ellipse after Ramsay (1967)

## Usage

``` r
mean_strain_ellipse_ramsay(
  r,
  phi = NULL,
  boot = TRUE,
  resamples = 1000,
  boot.values = FALSE
)
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

list. `R` aspect ratio of strain ellipse, `Ri` initial aspect ratio,
`Fl` the fluctuation angle (after Ramsay 1967), and `phi` the mean
orientation of the strain ellipse long axis. If `boot=TRUE`, then the
bootstrapped 95% confidence interval for the values are added. If
`boot.values=TRUE`, a matrix containing the bootstrapped R and phi
values are added.

## References

Ramsay (1976), Folding and Fracturing of Rocks, McGraw-Hill Book
Company.

Ramsay, J. G., & Huber, M. I. (1983). The Techniques of Modern
Structural Geology: Strain Analyses (Vol. 1). London: Academic Press.

## See also

[`mean_strain_ellipse()`](https://tobiste.github.io/structr/reference/mean_strain_ellipse.md)

## Examples

``` r
data(ramsay)
mean_strain_ellipse_ramsay(ramsay[, "R"], ramsay[, "phi"])
#> $R
#> [1] 1.850135
#> 
#> $Ri
#> [1] 1.762033
#> 
#> $Fl
#> [1] 65.79613
#> 
#> $phi
#> [1] 25.83207
#> 
#> $R_CI
#> [1] 1.633095 1.944634
#> 
#> $Ri_CI
#> [1] 1.540705 1.762033
#> 
#> $Fl_CI
#> [1] 46.67333 65.79613
#> 
#> $phi_CI
#> [1] 24.86458 26.85261
#> 
```
