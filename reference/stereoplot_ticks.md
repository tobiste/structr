# Stereoplot Tickmarks

Adds stereoplot rickmarks to an existing plot

## Usage

``` r
stereoplot_ticks(
  length = 0.02,
  angle = 10,
  labels = FALSE,
  ladj = 2 * length,
  radius = 1,
  rotation = 0,
  ...
)
```

## Arguments

- length:

  numeric. Length of ticks as fraction of `radius`

- angle:

  numeric. Division angle in degrees

- labels:

  this can either be a logical value specifying whether (numerical)
  annotations are to be made next to the tickmarks, or a character or
  expression vector of labels to be placed next to the tickpoints.

- ladj:

  adjustment for all labels away from origin of stereoplot circle. This
  essentially an amount that is added to `radius` and the length of the
  ticks.

- radius:

  numeric. Radius of circle

- rotation:

  numeric. Rotation (positive for counter-clockwise) of tickmarks and
  labels

- ...:

  optional arguments passed to
  [`graphics::segments()`](https://rdrr.io/r/graphics/segments.html) and
  [`graphics::text()`](https://rdrr.io/r/graphics/text.html)

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
[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
plot(c(-1, 1), c(-1, 1), type = "n", asp = 1)
stereoplot_frame()
stereoplot_ticks(length = 0.05, angle = 45, col = "blue", lwd = 2, labels = TRUE)
```
