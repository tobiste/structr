# Plot Bootstrapped Confidence Ellipse

Adds an ellipse marking the bootstrapped confidence interval of the
arithmetic mean to an existing plot

## Usage

``` r
stereo_confidence(
  x,
  params = list(),
  .center = TRUE,
  col = par("col"),
  cex = par("cex"),
  pch = 16,
  upper.hem = NULL,
  earea = NULL,
  radius = NULL,
  ...
)
```

## Arguments

- x:

  Spherical object or a list containing the output of an earlier call of
  [`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)

- params:

  list. Parameters passed to
  [`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)

- .center:

  logical. Whether the ellipse's center should be plotted?

- col:

  Color of the ellipse and its center

- pch, cex:

  Plotting symbol and size of the ellipse center. Ignored if `.center`
  is `FALSE`.

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

  graphical parameters passed to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html)

## Value

output of
[`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)

## See also

[`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)

Other stereo-plot:
[`arrows()`](https://tobiste.github.io/structr/reference/arrows.md),
[`fault-plot`](https://tobiste.github.io/structr/reference/fault-plot.md),
[`lines()`](https://tobiste.github.io/structr/reference/lines.md),
[`plot-spherical`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`points.spherical()`](https://tobiste.github.io/structr/reference/points.spherical.md),
[`stereo_cones`](https://tobiste.github.io/structr/reference/stereo_cones.md),
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
plot(example_lines, col = "grey")
stereo_confidence(example_lines, params = list(n = 100, res = 100), col = "red")
```
