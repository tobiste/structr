# Watson Distribution

A naive acceptance-rejection sampling algorithm, based on bounding the
density (with respect to the distance from `mu`) with a constant. For
large `kappa`, this method grows inefficient. For `kappa` == 100, about
13 tries are needed per success. For `kappa` == -100, about 18 tries are
needed.

## Usage

``` r
rwatson(n, mu, k)
```

## Source

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

## Value

vector of class `mu` of length `n`

## See also

Other random:
[`rbing`](https://tobiste.github.io/structr/reference/rbing.md),
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md),
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md),
[`rrot()`](https://tobiste.github.io/structr/reference/rrot.md),
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md),
[`vonmises-fisher`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)

## Examples

``` r
set.seed(20250411)
r <- rwatson(100, mu = Ray(120, 50), k = 10)

contour(r)
points(r)
```
