# Principal component (geodesic) analysis in the tangent space.

Principal component (geodesic) analysis in the tangent space.

## Usage

``` r
# S3 method for class 'Vec3'
prcomp(x, center = geodesic_mean(x), n = 0L)

# S3 method for class 'Ray'
prcomp(x, center = geodesic_mean(x), n = 0L)

# S3 method for class 'Line'
prcomp(x, center = geodesic_mean(x), n = 0L)

# S3 method for class 'Plane'
prcomp(x, center = geodesic_mean(x), n = 0L)

# S3 method for class 'Pair'
prcomp(x, center = geodesic_mean(x), n = 0L, group = NULL)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- center:

  A spherical object. Typically the geodesic mean of the x.

- n:

  real number (integer, 0 or \>= 3). The number of points to return on
  each of the four geodesics through the center.

- group:

  Symmetry group of `x`. See
  [`symmetry_group()`](https://tobiste.github.io/structr/reference/symmetry_group.md)
  for details If `NULL`, the group will be automatically picked based on
  the class of `x`.

## Value

A list consisting of

- `rotation`:

  3x3 rotation matrix

- `magnitudes`:

  2D real vector, non-negative. Magnitudes are analogous to sample
  standard deviations. They are in decreasing order and quantify how
  dispersed the data are in those two directions.

- `directions`:

  2x2 real matrix, whose columns are unit-length vectors. The
  corresponding directions to the magnitudes

- `pcsFromRay`:

  function to convert rays to 2-dimensional vectors

- `rayFromPCs`:

  function to convert 2-dimensional vectors to rays

- `curves`:

  list of two lists of (2 `n` + 1) rays, only if `n` \>= 1

- `tangents`:

  Tangents from `directions` and `rotation`

If `x` is a `"Pair"` object the list only contains `magnitudes`,
`directions` and `curves` (if `n`\>=1).

## Details

Appropriate only if the data set is tightly concentrated about the
center.

## Examples

``` r
res <- prcomp(example_lines, n = 10)

stereoplot(sub = paste("SD1:", round(res$magnitudes[1], 2), "| SD2:",round(res$magnitudes[2], 2)))
points(example_lines, pch = 16, cex = .8)
invisible(lapply(res$curves, stereo_lines, col = "red", lwd = 1.5))
```
