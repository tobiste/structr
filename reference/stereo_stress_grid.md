# Stress Grid

Plots the slip directions for a grid of fault plane normals and a given
stress tensor in a stereographic/equal-area projection

## Usage

``` r
stereo_stress_grid(sigma, n = 15, radius = NULL, ...)
```

## Arguments

- sigma:

  symmetric 3x3 matrix. The (reduced) stress tensor.

- n:

  integer. The size of the grid of fault planes. `20` by default.

- radius:

  numeric. Radius of circle. Defaults to `getOption("structr.radius")`.

- ...:

  optional arguments passed to
  [`stereo_hoeppener()`](https://tobiste.github.io/structr/reference/arrows.md)

## Value

object of class `"Fault"`

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md);
[`sigma2slip()`](https://tobiste.github.io/structr/reference/sigma2slip.md)

## Examples

``` r
sigma <- reduced_stress(angelier1990$AVB)

dev.off()
#> null device 
#>           1 
stereoplot(guides = FALSE)
stereo_stress_grid(sigma, type = 'hoeppener')
#> Warning: graphical parameter "type" is obsolete

pr <- sigma2stress(sigma)
points(pr$axes, col = 2:4, pch = 16)
text(pr$axes, label = rownames(pr$axes), col = 2:4, adj = -.25)
```
