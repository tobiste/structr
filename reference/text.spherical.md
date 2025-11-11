# Add Points to a Plot

Add Points to a Plot

## Usage

``` r
# S3 method for class 'spherical'
text(x, labels = seq_along(x[, 1]), upper.hem = FALSE, earea = TRUE, ...)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- labels:

  a character vector or
  [expression](https://rdrr.io/r/base/expression.html) specifying the
  *text* to be written. An attempt is made to coerce other language
  objects (names and calls) to expressions, and vectors and other
  classed objects to character vectors by
  [`as.character`](https://rdrr.io/r/base/character.html). If `labels`
  is longer than `x` and `y`, the coordinates are recycled to the length
  of `labels`.

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`, the default).

- earea:

  logical `TRUE` for Lambert equal-area projection (also "Schmidt net";
  the default), or `FALSE` for meridional stereographic projection (also
  "Wulff net" or "Stereonet").

- ...:

  arguments passed to
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
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md)

## Examples

``` r
stereoplot()
points(Line(c(90, 80), c(10, 75)), col = 1:2)
text(Line(c(90, 80), c(10, 75)), labels = c("L1", "L2"), col = 1:2, pos = 3)
```
