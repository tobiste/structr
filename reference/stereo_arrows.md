# Add Arrows to a Stereoplot

A quiver plot displays displacement vectors into pointing into the
direction of movement.

## Usage

``` r
stereo_arrows(
  x,
  sense,
  scale = 0.1,
  angle = 10,
  length = 0.1,
  upper.hem = FALSE,
  earea = TRUE,
  ...
)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- sense:

  numeric. Sense of the line on a fault plane. Either `1`or `-1` for
  normal or thrust offset, respectively. The "sense" is the sign of the
  fault's rake (see
  [`Fault_from_rake()`](https://tobiste.github.io/structr/reference/fault_from_rake.md)
  for details).

- scale:

  numeric. Scales the length of the vector. `0.1` by default

- angle:

  numeric. Angle from the shaft of the arrow to the edge of the arrow
  head.

- length:

  numeric. Length of the edges of the arrow head (in inches).

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`, the default).

- earea:

  logical `TRUE` for Lambert equal-area projection (also "Schmidt net";
  the default), or `FALSE` for meridional stereographic projection (also
  "Wulff net" or "Stereonet").

- ...:

  arguments passed to
  [`graphics::arrows()`](https://rdrr.io/r/graphics/arrows.html)

## See also

[`hoeppener()`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`angelier()`](https://tobiste.github.io/structr/reference/fault-plot.md)

Other stereo-plot:
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines.spherical()`](https://tobiste.github.io/structr/reference/lines.spherical.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
[`stereo_cones`](https://tobiste.github.io/structr/reference/stereo_cones.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md),
[`stereo_contour`](https://tobiste.github.io/structr/reference/stereo_contour.md),
[`stereo_lines()`](https://tobiste.github.io/structr/reference/stereo_lines.md),
[`stereo_point()`](https://tobiste.github.io/structr/reference/stereo_point.md),
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
set.seed(20250411)
stereoplot()
p <- rvmf(n = 100)
points(p, pch = 16, cex = .5)
stereo_arrows(p, sense = 1, col = "red")
```
