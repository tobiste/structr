# von Mises-Fisher Distribution

Density and random generation for the spherical normal distribution with
mean and concentration parameter (\\\kappa\\) .

## Usage

``` r
rvmf(
  n = 100,
  mu = Vec3(1, 0, 0),
  k = 5,
  method = c("geologyGeometry", "rotasym")
)

dvmf(x, mu, k = 5)
```

## Source

Adapted fom
[`rotasym::r_vMF()`](https://rdrr.io/pkg/rotasym/man/vMF.html) and
[`rotasym::d_vMF()`](https://rdrr.io/pkg/rotasym/man/vMF.html), and
`geologyGeometry` by Davis, J.R.

## Arguments

- n:

  integer. number of random samples to be generated

- mu:

  Mean vector. object of class `"Vec3"`, `"Line"`, `"Ray"`, or
  `"Plane"`, where the rows are the observations and the columns are the
  coordinates.

- k:

  numeric. The concentration parameter (\\\kappa\\) of the von
  Mises-Fisher distribution

- method:

  character. Algorithm to generate random vectors from a Fisher
  distribution. Either `"geologyGeometry"` (the default) to pick the
  `rayFisher()` algorithm from the *geologyGeometry* code compilation,
  or `"rotasym"` to pick the
  [`rotasym::r_vMF()`](https://rdrr.io/pkg/rotasym/man/vMF.html)
  algorithm from the *rotasym* package.

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

## See also

Other random:
[`rbing`](https://tobiste.github.io/structr/reference/rbing.md),
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md),
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md),
[`rrot()`](https://tobiste.github.io/structr/reference/rrot.md),
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md),
[`rwatson()`](https://tobiste.github.io/structr/reference/rwatson.md)

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Ray(120, 50), k = 5)
contour(x)
points(x)


dx <- dvmf(x, mu = Ray(120, 50))
head(dx)
#>            [,1]
#> [1,] 0.65684755
#> [2,] 0.14268609
#> [3,] 0.34589805
#> [4,] 0.05253657
#> [5,] 0.67050140
#> [6,] 0.18688448

plot(x, col = assign_col(dx))
```
