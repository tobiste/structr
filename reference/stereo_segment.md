# Great-circle Segment Between Two Vectors

Plots the great-circle segment between two vectors

## Usage

``` r
stereo_segment(
  x,
  y,
  upper.hem = NULL,
  earea = NULL,
  n = 100L,
  radius = NULL,
  ...
)
```

## Arguments

- x, y:

  objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`). Defaults to
  `getOption("structr.upper.hem")`.

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection,
  or `FALSE` for meridional stereographic projection. Defaults to
  `getOption("structr.earea")`.

- n:

  integer. number of points along great-circle (100 by default)

- radius:

  numeric. Radius of circle. Defaults to `getOption("structr.radius")`.

- ...:

  optional graphical parameters passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## See also

[`slerp()`](https://tobiste.github.io/structr/reference/slerp.md),
[stereo_greatcircle](https://tobiste.github.io/structr/reference/stereo_cones.md),
[`stereo_lines()`](https://tobiste.github.io/structr/reference/stereo_lines.md)

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
[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`stereoplot_guides()`](https://tobiste.github.io/structr/reference/stereoplot_guides.md),
[`stereoplot_ticks()`](https://tobiste.github.io/structr/reference/stereoplot_ticks.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md)

## Examples

``` r
x <- Line(120, 7)
y <- Line(10, 13)
plot(rbind(x, y))
stereo_segment(x, y, col = "red")


# For multiple segments use lapply():
set.seed(20250411)
mu <- Line(45, 10)
x <- rvmf(100, mu = mu) |> Line()
plot(x)
invisible(lapply(seq_len(nrow(x)), FUN = function(i) {
  stereo_segment(x[i, ], mu, col = i)
}))
points(mu, pch = 16, col = "white")
```
