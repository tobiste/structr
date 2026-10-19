# Spherical Projection of Cones

Visualization of small-circles and great-circles in a stereographic or
equal-area projection.

## Usage

``` r
stereo_smallcircle(
  x,
  d = 90,
  col = par("col"),
  N = 1000,
  upper.hem = NULL,
  earea = NULL,
  lty = par("lty"),
  lwd = par("lwd"),
  fill = FALSE,
  border = NA,
  radius = NULL,
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
  or lower hemisphere (`FALSE`). Defaults to
  `getOption("structr.upper.hem")`.

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection,
  or `FALSE` for meridional stereographic projection. Defaults to
  `getOption("structr.earea")`.

- fill:

  logical. Whether to fill the inner part of the small-circle? `FALSE`
  by default.

- border:

  Color of the filled small-circle's outline (ignored if `fill=FALSE`)

- radius:

  numeric. Radius of circle. Defaults to `getOption("structr.radius")`.

- ...:

  optional graphical parameters passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html) and (if
  `fill=TRUE`)
  [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html)

## See also

[`lines.spherical()`](https://tobiste.github.io/structr/reference/lines.md),
[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`stereo_lines()`](https://tobiste.github.io/structr/reference/stereo_lines.md)

Other stereo-plot:
[`arrows()`](https://tobiste.github.io/structr/reference/arrows.md),
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines()`](https://tobiste.github.io/structr/reference/lines.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
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
stereo_point(Line(90, 5), lab = "L")
stereo_smallcircle(Line(90, 5), d = 10)
stereo_point(Plane(120, 30), lab = "P", col = "red")
stereo_greatcircle(Plane(120, 30), col = "red")


stereoplot()
stereo_point(Line(c(129, 90), c(30, 5)), lab = c("L1", "L2"))
stereo_smallcircle(Line(c(129, 90), c(30, 5)), d = c(10, 5), 
  col = 1:2, lty = 1:2, lwd = 1:2)

#> [[1]]
#> NULL
#> 
#> [[2]]
#> NULL
#> 

# Filled cones:
stereoplot()
stereo_smallcircle(Line(c(90, 120), c(5, 5)), d = c(5, 20), 
  col = c('grey60', 'grey40'), border = c('red', 'blue'), fill = TRUE)
stereo_point(Line(c(90, 120), c(5, 5)), col = c('red', 'blue'))
```
