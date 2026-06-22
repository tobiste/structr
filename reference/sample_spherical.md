# Random Samples and Permutations of Spherical Objects

`sample_spherical` takes a sample of the specified size from the
elements of the spherical object `x` using either with or without
replacement.

## Usage

``` r
sample_spherical(x, size, replace = FALSE, prob = NULL)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`.

- size:

  a non-negative integer giving the number of items to choose.

- replace:

  should sampling be with replacement?

- prob:

  a vector of probability weights for obtaining the elements of the
  vector being sampled.

## Value

A spherical object of length `size` with elements drawn from `x`

## Examples

``` r
set.seed(20250411)
x <- rvmf(n = 100, mu = Line(90, 45))
sample_spherical(x, size = 5)
#> Ray object (n = 5):
#>        azimuth    plunge
#> [1,] 124.88841 74.965893
#> [2,]  58.03926 24.560547
#> [3,]  90.03802  5.859954
#> [4,]  58.58939  9.788074
#> [5,] 121.90908 17.538452
```
