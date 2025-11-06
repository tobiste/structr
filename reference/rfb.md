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

[`rvmf()`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)
to draw samples from the von Mises Fisher distribution around a
specified mean vector.
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md) to
draw from a Kent-distribution.
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md)
to draw from a a spherical uniform distribution.

## Examples

``` r
set.seed(20250411)
x <- rfb(100, mu = Line(120, 50), k = 5, A = diag(c(-1, 0, 1)))
plot(x)
```
