# Stereoplot frame

Adds a (primitive) circle with given radius to an existing plot

## Usage

``` r
stereoplot_frame(n = 512L, radius = 1, ...)
```

## Arguments

- n:

  integer. Resolution of circle's line

- radius:

  numeric. Radius of circle

- ...:

  optional arguments passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## See also

[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md)

## Examples

``` r
plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
stereoplot_frame(col = "red", lwd = 3)
```
