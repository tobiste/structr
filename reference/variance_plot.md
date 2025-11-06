# Variance Plot

Shows the greatcircle of the shortest distance between a set of vectors
to a specified vector in a stereoplot. The greatcircles are color-coded
by the angular distance.

## Usage

``` r
variance_plot(
  x,
  y = NULL,
  .mean = c("geodesic", "arithmetic", "projected"),
  show.center = TRUE,
  segments = TRUE,
  upper.hem = FALSE,
  earea = TRUE,
  ...
)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- y:

  The vector from which the variance should be visualized (only one
  vector allowed). When `NULL`, then the mean vector of `x` is used (the
  default).

- .mean:

  character. The type of mean to be used if `y` is `NULL`. One of
  `"geodesic"` (the default), `"arithmetic"` or `"projected"`.

- show.center:

  logical. Whether the center point `y` of the mean should be
  highlighted in the plot.

- segments:

  logical. Whether the segments should be shown or only the points?

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`, the default).

- earea:

  logical `TRUE` for Lambert equal-area projection (also "Schmidt net";
  the default), or `FALSE` for meridional stereographic projection (also
  "Wulff net" or "Stereonet").

- ...:

  optional arguments passed to
  [`assign_col()`](https://tobiste.github.io/structr/reference/colorize.md)

## Value

list. `angles` is a vector of the geodesic angles (in degrees) between
all vectors in `x` and `y` (or the mean), and `var` is a scalar giving
the Fréchet variance.

## See also

[`stereo_segment()`](https://tobiste.github.io/structr/reference/stereo_segment.md),
[`sph_mean()`](https://tobiste.github.io/structr/reference/stats.md),
[`geodesic_mean()`](https://tobiste.github.io/structr/reference/geodesic-mean.md),
[`projected_mean()`](https://tobiste.github.io/structr/reference/projected_mean.md),
[`geodesic_var()`](https://tobiste.github.io/structr/reference/geodesic-var.md)

## Examples

``` r
variance_plot(example_lines)

variance_plot(example_planes, example_planes[1, ], segments = FALSE)
```
