# Stereographic Projection

Initialize the plot for equal-area stereographic projections (Wulff) or
Lambert Equal-Area projections (Schmidt).

## Usage

``` r
stereoplot(
  earea = TRUE,
  guides = TRUE,
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
  radius = 1
)
```

## Source

Adapted from the `RFOC` package

## Arguments

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection
  (the default), or `FALSE` for meridional stereographic projection.

- guides:

  logical. Whether guides should be added to the plot (`TRUE` by
  default)

- d:

  integer. Angle distance between guides. Default: 10

- col:

  color of guide lines

- lwd:

  linewidth of guide lines

- lty:

  linetype of guide lines

- border.col:

  color of primitive circle (frame), center-cross and ticks of the
  stereo plot

- title, sub:

  character. Title and subtitle of plot

- origin.text:

  character. Text at origin of stereoplot.

- labels:

  this can either be a logical value specifying whether (numerical)
  annotations are to be made next to the tickmarks, or a character or
  expression vector of labels to be placed next to the tickpoints.

- ladj:

  adjustment for all labels away from origin of stereoplot circle. This
  essentially an amount that is added to `radius` and the length of the
  ticks.

- centercross:

  logical. Whether a center cross should be added (`TRUE` by default)

- ticks:

  integer. Angle between ticks. if `NULL` (the default), no ticks are
  drawn.

- radius:

  numeric. Radius of circle

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
[`stereo_point()`](https://tobiste.github.io/structr/reference/stereo_point.md),
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
stereoplot()


stereoplot(ticks = 30, title = "title", sub = "subtitle", border.col = "purple", labels = TRUE)
```
