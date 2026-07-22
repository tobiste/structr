# Center gridlines on a given point

Center gridlines on a given point

## Usage

``` r
rotate_stereogrid(x, d = 10, col = "gray90", lwd = 0.5, lty = 1, ...)
```

## Arguments

- x:

  center position of grid lines.

- d:

  integer. Angle distance between guides. Defaults to
  `getOption("structr.d")`.

- col:

  Color of guide lines. Defaults to `getOption("structr.col")`.

- lwd:

  Width of guide lines. Defaults to `getOption("structr.lwd")`.

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
