# Lines in a Stereoplot

Draws simple lines between vector points in stereographic or equal-area
projection

## Usage

``` r
stereo_lines(x, upper.hem = NULL, earea = NULL, radius = NULL, ...)
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

- radius:

  numeric. Radius of circle. Defaults to `getOption("structr.radius")`.

- ...:

  optional graphical parameters passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## Value

two-column matrix of the stereographic or equal-area coordinates

## See also

[`slerp()`](https://tobiste.github.io/structr/reference/slerp.md),
[stereo_greatcircle](https://tobiste.github.io/structr/reference/stereo_cones.md),
`stereo_lines()`,
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md)

Other stereo-plot:
[`arrows()`](https://tobiste.github.io/structr/reference/arrows.md),
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines()`](https://tobiste.github.io/structr/reference/lines.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
[`stereo_cones`](https://tobiste.github.io/structr/reference/stereo_cones.md),
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
plot(example_lines, col = "grey")
stereo_lines(example_lines[1:2, ], col = "red")
```
