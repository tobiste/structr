# Flinn diagram

Flinn diagram

## Usage

``` r
flinn_plot(
  x,
  main = "Flinn diagram",
  R.max = NULL,
  log = FALSE,
  add = FALSE,
  ...
)

# Default S3 method
flinn_plot(
  x,
  main = "Flinn diagram",
  R.max = NULL,
  log = FALSE,
  add = FALSE,
  ...
)

# S3 method for class 'ortensor'
flinn_plot(x, ...)

# S3 method for class 'ellipsoid'
flinn_plot(x, ...)

# S3 method for class 'spherical'
flinn_plot(x, ...)
```

## Arguments

- x:

  accepts the following objects: a two-column matrix where first column
  is the ratio of maximum strain and intermediate strain (X/Y) and
  second column is the the ratio of intermediate strain and minimum
  strain (Y/Z); objects of class `"Vec3"`, `"Line"`, `"Ray"`, or
  `"Plane"`; or `"ortensor"` objects.

- main:

  character. The main title (on top).

- R.max:

  numeric. Maximum aspect ratio for scaling.

- log:

  logical. Whether the axes should be in logarithmic scale.

- add:

  logical. Should data be plotted to an existing plot?

- ...:

  plotting arguments passed to
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html)

## Value

plot and when stored as an object, the multiplication factors for X, Y
and Z.

## References

Flinn, D. (1965). On the Symmetry Principle and the Deformation
Ellipsoid. Geological Magazine, 102(1), 36–45.
[doi:10.1017/S0016756800053851](https://doi.org/10.1017/S0016756800053851)

## See also

Other fabric-plot:
[`hsu_plot()`](https://tobiste.github.io/structr/reference/hsu_plot.md),
[`vollmer-plot`](https://tobiste.github.io/structr/reference/vollmer-plot.md),
[`woodcock_plot()`](https://tobiste.github.io/structr/reference/woodcock_plot.md)

## Examples

``` r
data(holst)
R_XY <- holst[, "R_XY"]
R_YZ <- holst[, "R_YZ"]
flinn_plot(cbind(R_XY, R_YZ), log = FALSE, col = "#B63679", pch = 16)

flinn_plot(cbind(R_XY, R_YZ), log = TRUE, col = "#B63679", pch = 16, type = "b")


set.seed(20250411)
mu <- Line(120, 50)
x <- rvmf(100, mu = mu, k = 1)
flinn_plot(x, R.max = 2)

set.seed(20250411)
y <- rvmf(100, mu = mu, k = 20)
flinn_plot(ortensor(y), col = "red", R.max = 2, add = TRUE)
```
