# Spherical Fisher-Bingham distribution

Simulation of random values from a spherical Fisher-Bingham
distribution.

## Usage

``` r
rfb(n = 100, mu = Vec3(1, 0, 0), k = 5, A)
```

## Source

Adapted from
[`Directional::rfb()`](https://rdrr.io/pkg/Directional/man/rfb.html)

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

- A:

  symmetric matrix

## See also

Other random:
[`rbing`](https://tobiste.github.io/structr/reference/rbing.md),
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md),
[`rrot()`](https://tobiste.github.io/structr/reference/rrot.md),
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md),
[`vonmises-fisher`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)

## Examples

``` r
set.seed(20250411)
x <- rfb(100, mu = Line(120, 50), k = 5, A = diag(c(-1, 0, 1)))
contour(x)
```
