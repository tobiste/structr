# Mohr Circle Diagram for Strain

Plots the Mohr Circle for Strain

## Usage

``` r
Mohr_strain(
  lambda1,
  lambda2 = NA,
  lambda3,
  phi = NULL,
  fg = par("col"),
  bg = "lightgray",
  fg23 = par("col"),
  bg23 = "white",
  fg12 = par("col"),
  bg12 = "white",
  axes = TRUE,
  col = "black",
  n = 512,
  full.circle = FALSE,
  include.zero = TRUE,
  xlim = NULL,
  ylim = NULL,
  round = 1,
  ...
)
```

## Arguments

- lambda1, lambda2, lambda3:

  numeric. Magnitude of quadratic elongation \\\lambda\\

- phi:

  numeric. (optional) Angle (in degrees) for a specific strain

- fg, fg12, fg23:

  border color for the Mohr Circles spanning lambda1-lambda3,
  lambda1-lambda2, and lambda2-lambda3, respectively

- bg, bg12, bg23:

  fill color for the Mohr Circles spanning lambda1-lambda3,
  lambda1-lambda2, and lambda2-lambda3, respectively

- axes:

  logical. Show axis?

- col:

  color for the stress state for a given `phi`.

- n:

  integer. Resolution given amount of points along the generated path
  representing the full Mohr circle (`512` by default).

- full.circle:

  logical. Should the complete Mohr circle be shown, or only the upper
  (positive shear stress) part of the circle?

- include.zero:

  logical. the plot range be extended to include `lambda = 0`?

- xlim, ylim:

  range of plot

- round:

  integer indicating the number of decimal places to be used for
  rounding.

- ...:

  optional graphical parameters.

## See also

[`Mohr_plot()`](https://tobiste.github.io/structr/reference/Mohr_plot.md)
for Stress.
[strain](https://tobiste.github.io/structr/reference/strain.md) for
converting strain quantities

## Examples

``` r
Mohr_strain(lambda1 = 4, lambda3 = 0.25, phi = 25, col = 'red')

Mohr_strain(lambda1 = 4, lambda2 = 1, lambda3 = 0.25, phi = 25, col = 'red', full.circle = TRUE, axes = FALSE)
```
