# Center grid lines on a given vector

Center grid lines on a given vector

## Usage

``` r
rotate_stereogrid(
  x,
  d = 10,
  col = "gray90",
  lwd = 0.5,
  equator = TRUE,
  equator_lwd = 1.5 * lwd,
  lty = 1,
  ...
)
```

## Arguments

- x:

  center position of grid lines.

- d:

  angle spacing between grid lines of the projection

- col:

  Color of guide lines. Defaults to `getOption("structr.col")`.

- lwd:

  Width of guide lines. Defaults to `getOption("structr.lwd")`.

- equator:

  logical. Whether the grid equator should be shown no matter how the
  grid is constructed via `d`. `TRUE` by default.

- equator_lwd:

  numeric. The line width of the drawn equator grid line. By default,
  the equator is `1.5` times thicker then `lwd`.

- lty:

  Type of guide lines. Defaults to `getOption("structr.lty")`.

- ...:

  arguments passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## Examples

``` r
stereoplot(guide = FALSE)
rotate_stereogrid(Plane(120, 50), earea = FALSE)
```
