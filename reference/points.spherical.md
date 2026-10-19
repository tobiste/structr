# Add Points to a Plot

Add Points to a Plot

## Usage

``` r
# S3 method for class 'spherical'
points(x, upper.hem = NULL, earea = NULL, ...)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`). Defaults to
  `getOption("structr.upper.hem")`.

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection,
  or `FALSE` for meridional stereographic projection. Defaults to
  `getOption("structr.earea")`.

- ...:

  arguments passed to
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html)

## See also

Other stereo-plot:
[`arrows()`](https://tobiste.github.io/structr/reference/arrows.md),
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines()`](https://tobiste.github.io/structr/reference/lines.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
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
stereoplot()
points(rvmf(n = 100))
points(Plane(120, 30), col = "red", pch = 19)
```
