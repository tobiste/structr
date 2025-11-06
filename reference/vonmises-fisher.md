# von Mises-Fisher Distribution

Density and random generation for the spherical normal distribution with
mean and concentration parameter (\\\kappa\\) .

## Usage

``` r
rvmf(n = 100, mu = Vec3(1, 0, 0), k = 5)

dvmf(x, mu, k = 5)
```

## Source

Adapted fom
[`rotasym::r_vMF()`](https://rdrr.io/pkg/rotasym/man/vMF.html) and
[`rotasym::d_vMF()`](https://rdrr.io/pkg/rotasym/man/vMF.html)

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

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

## See also

[`runif.spherical()`](https://tobiste.github.io/structr/reference/runif.spherical.md)
for alternative algorithms to generate uniform distributed samples on a
sphere,
[`rkent()`](https://tobiste.github.io/structr/reference/rkent.md) for
Kent distribution,
[`rfb()`](https://tobiste.github.io/structr/reference/rfb.md) for
Fisher-Bingham distribution.

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Line(120, 50), k = 5)
dx <- dvmf(x, mu = Line(120, 50))
head(dx)
#>           [,1]
#> [1,] 0.3539639
#> [2,] 0.3456839
#> [3,] 0.3298219
#> [4,] 0.2878725
#> [5,] 0.3416724
#> [6,] 0.3087497

plot(x, col = assign_col(dx))
```
