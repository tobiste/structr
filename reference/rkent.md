# Kent Distribution

Simulation of random values from a spherical Kent distribution.

## Usage

``` r
rkent(n, mu = Vec3(1, 0, 0), k = 5, b)
```

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

- b:

  numeric. \\\beta\\ (ellipticity): \\0 \leq \beta \< \kappa\\

## See also

Other random:
[`rbing`](https://tobiste.github.io/structr/reference/rbing.md),
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md),
[`rrot()`](https://tobiste.github.io/structr/reference/rrot.md),
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md),
[`rwatson()`](https://tobiste.github.io/structr/reference/rwatson.md),
[`vonmises-fisher`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)

## Examples

``` r
set.seed(20250411)
r <- rkent(100, mu = Ray(120, 50), k = 5, b = 1)

contour(r)
points(r)
```
