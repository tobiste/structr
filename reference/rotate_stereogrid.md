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

  integer. Angle distance between guides. Default: 10

- col:

  Color of guide lines

- lwd:

  Width of guide lines

- lty:

  Type of guide lines

- ...:

  arguments passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## Examples

``` r
stereoplot(guide = FALSE)
rotate_stereogrid(Plane(120, 50), earea = FALSE)
```
