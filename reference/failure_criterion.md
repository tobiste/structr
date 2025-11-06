# Failure Criterion

Adds the Griffith-Coulomb-fracture criterion to a plot.

## Usage

``` r
failure_criterion(
  sigma_n = seq(-100, 100, 0.01),
  cohesion = 1,
  friction = 0.6,
  ...
)
```

## Arguments

- sigma_n:

  numeric. A vector of normal stresses, for which the critical shear
  stress should be calculated.

- cohesion:

  numeric. Cohesion

- friction:

  numeric. Coefficient of friction

- ...:

  optional plotting arguments passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## Value

a matrix of the normal and shear stresses.

## Examples

``` r
Mohr_plot(sigma1 = 66, sigma2 = 30, sigma3 = 20, xlim = c(-50, 125), full.circle = TRUE)
failure_criterion(col = "red")
```
