# Stereographic and Equal Area Projection

Initialize the plot for equal-angle stereographic projections (Wulff) or
Lambert equal-area projections (Schmidt).

## Usage

``` r
stereoplot(
  earea = NULL,
  guides = NULL,
  d = 10,
  col = "gray90",
  lwd = 0.5,
  lty = 1,
  border.col = "black",
  title = NULL,
  sub = NULL,
  origin.text = "N",
  labels = FALSE,
  ladj = 0.05,
  centercross = TRUE,
  ticks = 90,
  radius = NULL,
  center = NULL
)
```

## Source

Adapted from the `RFOC` package

## Arguments

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection,
  or `FALSE` for meridional stereographic projection. Defaults to
  `getOption("structr.earea")`.

- guides:

  logical. Whether guides should be added to the plot. Defaults to
  `getOption("structr.guides")`.

- d:

  integer. Angle distance between guides. Default: 10

- col:

  Color of guide lines

- lwd:

  Width of guide lines

- lty:

  Type of guide lines

- border.col:

  color of primitive circle (frame), center-cross and ticks of the
  stereo plot

- title, sub:

  character. Title and subtitle of plot

- origin.text:

  character. Text at origin of plot

- labels:

  this can either be a logical value specifying whether (numerical)
  annotations are to be made next to the tick marks, or a character or
  expression vector of labels to be placed next to the tick points.

- ladj:

  adjustment for all labels away from origin of projection circle. This
  essentially an amount that is added to `radius` and the length of the
  ticks.

- centercross:

  logical. Whether a center cross should be added (`TRUE` by default)

- ticks:

  integer. Angle between ticks. if `NULL` (the default), no ticks are
  drawn.

- radius:

  numeric. Radius of circle. Defaults to `getOption("structr.radius")`.

- center:

  An object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`
  specifying the center of the projection If `NULL` (the default), the
  center is at the origin of the plot.

## See also

[structr-options](https://tobiste.github.io/structr/reference/structr-options.md)

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
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
stereoplot()


stereoplot(ticks = 30, title = "title", sub = "subtitle", border.col = "purple", labels = TRUE)


stereoplot(center = Line(120, 50))
```
