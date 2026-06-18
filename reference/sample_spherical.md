# Random Samples and Permutations of Spherical Objects

Random Samples and Permutations of Spherical Objects

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

## Examples

``` r
set.seed(20250411)
x <- rvmf(n = 100, mu = Line(90, 45))
sample_spherical(x, size = 5)
#> Line object (n = 5):
#>        azimuth   plunge
#> [1,]  87.14959 36.28062
#> [2,] 169.76083 18.63444
#> [3,]  82.39226 43.89736
#> [4,] 111.35962 45.75890
#> [5,]  85.78076 10.87919
```
