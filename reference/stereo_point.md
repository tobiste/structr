# Stereographic Projection of Lines and Planes

Visualization of lines, planes in a stereographic projection.

## Usage

``` r
stereo_point(
  x,
  col = 1,
  pch = 20,
  lab = NULL,
  text.pos = 4,
  cex = 1,
  upper.hem = FALSE,
  earea = TRUE,
  ...
)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- col:

  color

- pch:

  plotting character

- lab:

  character. text labels

- text.pos:

  position for labels

- cex:

  character expansion of labels

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`, the default).

- earea:

  logical `TRUE` for Lambert equal-area projection (also "Schmidt net";
  the default), or `FALSE` for meridional stereographic projection (also
  "Wulff net" or "Stereonet").

- ...:

  optional graphical parameters

## Value

two-column matrix of the stereographic or equal-area coordinates

## Note

`"Plane"` and `"Fault"` objects will be displayed as pole to the plane
(only).

## See also

Other stereo-plot:
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines.spherical()`](https://tobiste.github.io/structr/reference/lines.spherical.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
[`stereo_arrows()`](https://tobiste.github.io/structr/reference/stereo_arrows.md),
[`stereo_cones`](https://tobiste.github.io/structr/reference/stereo_cones.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md),
[`stereo_contour`](https://tobiste.github.io/structr/reference/stereo_contour.md),
[`stereo_lines()`](https://tobiste.github.io/structr/reference/stereo_lines.md),
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
stereoplot()
stereo_point(Line(c(90, 80), c(10, 75)), lab = c("L1", "L2"))
stereo_point(Plane(120, 30), lab = "P", col = "red")
```
