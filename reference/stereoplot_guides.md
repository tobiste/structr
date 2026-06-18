# Stereoplot Gridlines

Adds equal-area or equal-angle projection gridlines to an existing
stereoplot.

## Usage

``` r
stereoplot_guides(d = 10, earea = TRUE, radius = 1, center = NULL, ...)
```

## Arguments

- d:

  angle between grid lines

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection
  (the default), or `FALSE` for meridional stereographic projection.

- radius:

  numeric. Radius of circle

- center:

  Center position of gridlines. If `NULL` (the default), gridlines are
  centered on the origin of the stereoplot. Otherwise, this should is an
  spherical object.

- ...:

  optional arguments passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## See also

Other stereo-plot:
[`arrows()`](https://tobiste.github.io/structr/reference/arrows.md),
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines()`](https://tobiste.github.io/structr/reference/lines.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
[`stereo_cones`](https://tobiste.github.io/structr/reference/stereo_cones.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md),
[`stereo_contour`](https://tobiste.github.io/structr/reference/stereo_contour.md),
[`stereo_lines()`](https://tobiste.github.io/structr/reference/stereo_lines.md),
[`stereo_point()`](https://tobiste.github.io/structr/reference/stereo_point.md),
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
stereoplot_guides(d = 5, earea = FALSE, col = "green", rotation = 20)


plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
stereoplot_guides(d = 15, earea = TRUE, col = "orange", rotation = 90)


plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
stereoplot_guides(d = 15, earea = FALSE, center = Line(120, 50))
```
