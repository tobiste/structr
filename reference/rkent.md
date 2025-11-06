# Kent distribution

Simulation of random values from a spherical Kent distribution.

## Usage

``` r
rkent(n = 100, mu = Vec3(1, 0, 0), k = 5, b)
```

## Source

Adapted from
[`Directional::rkent()`](https://rdrr.io/pkg/Directional/man/rkent.html)

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

[`rvmf()`](https://tobiste.github.io/structr/reference/vonmises-fisher.md)
to draw samples from the von Mises Fisher distribution around a
specified mean vector.
[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md)
to draw from a a spherical uniform distribution.
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md) to draw
from a Fisher-Bingham distribution.

## Examples

``` r
set.seed(20250411)
x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
plot(x)
```
