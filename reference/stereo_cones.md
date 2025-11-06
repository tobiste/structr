# Stereographic Projection of Cones

Visualization of smallcircles and greatcircles in a stereographic
projection.

## Usage

``` r
stereo_smallcircle(
  x,
  d = 90,
  col = 1,
  N = 1000,
  upper.hem = FALSE,
  earea = TRUE,
  lty = 1,
  lwd = 1,
  BALL.radius = 1,
  ...
)

stereo_greatcircle(x, ...)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- d:

  numeric. conical angle in degrees.

- col, lty, lwd:

  color, line type, and line width parameters

- N:

  integer. number of points to calculate

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`, the default).

- earea:

  logical `TRUE` for Lambert equal-area projection (also "Schmidt net";
  the default), or `FALSE` for meridional stereographic projection (also
  "Wulff net" or "Stereonet").

- BALL.radius:

  numeric size of sphere

- ...:

  optional graphical parameters

## See also

Other stereo-plot:
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines.spherical()`](https://tobiste.github.io/structr/reference/lines.spherical.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
[`stereo_arrows()`](https://tobiste.github.io/structr/reference/stereo_arrows.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md),
[`stereo_contour`](https://tobiste.github.io/structr/reference/stereo_contour.md),
[`stereo_point()`](https://tobiste.github.io/structr/reference/stereo_point.md),
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
stereoplot()
stereo_point(Line(90, 5), lab = "L")
stereo_smallcircle(Line(90, 5), d = 10)
stereo_point(Plane(120, 30), lab = "P", col = "red")
stereo_greatcircle(Plane(120, 30), col = "red")


stereoplot()
stereo_point(Line(c(129, 90), c(30, 5)), lab = c("L1", "L2"))
stereo_smallcircle(Line(c(129, 90), c(30, 5)), d = c(10, 5), col = 1:2, lty = 1:2, lwd = 1:2)
```
