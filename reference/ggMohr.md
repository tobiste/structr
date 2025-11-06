# Mohr Circle plot

Plots the Mohr Circle diagram using ggplot functionality.

## Usage

``` r
ggMohr(
  sigma1,
  sigma2,
  sigma3 = NULL,
  units = "MPa",
  coulomb = c(1, 0.6),
  sliding = 0.81,
  fill = "gray",
  alpha = 0.5,
  show.info = TRUE,
  ...
)
```

## Arguments

- sigma1, sigma2, sigma3:

  numeric. Magnitudes of major, intermediate, and minor principal
  stresses. If only two principal stresses are given, only one Mohr
  Circle will be drawn, otherwise three.

- units:

  units of `sigma1`, `sigma2`, `sigma3` and cohesion (`"MPa"` by
  default). `NULL` if unitless (e.g. deviatoric stresses only)

- coulomb:

  numeric 2 element vector. Coulomb criterion containing the cohesion
  and the coefficient of sliding: (`c(70, 0.6)`)

- sliding:

  Sliding criteria (`0.81` by default)

- fill:

  fill color of Mohr circle spanned by sigma1 and sigma3

- alpha:

  opacity of Mohr circle spanned sigma1 and sigma3

- show.info:

  logical. Should Mohr parameters (mean stress and differential stress)
  be shown in the plot's caption?

- ...:

  optional parameters passed to
  [`ggforce::geom_circle()`](https://ggforce.data-imaginist.com/reference/geom_circle.html)

## Details

The subcaption gives the mean stress \\\sigma_m\\, the differential
stress \\\sigma_d\\, and the fracture angles \\\theta_f\\ and \\\alpha_f
= 90 - \theta_f\\

## See also

[`Mohr_plot()`](https://tobiste.github.io/structr/reference/Mohr_plot.md)

## Examples

``` r
ggMohr(1025, 450, 250)


# Unitless Mohr circle
ggMohr(5, 2, 0, units = NULL)
```
