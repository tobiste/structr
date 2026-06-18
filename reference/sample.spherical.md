# Random Samples and Permutations of Spherical Objects

Random Samples and Permutations of Spherical Objects

## Usage

``` r
# S3 method for class 'spherical'
sample(x, size, replace = FALSE, prob = NULL)
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

object of class `x`

## Examples

``` r
set.seed(20250411)
x <- rvmf(n = 100, mu = Line(90, 45))
sample(x, size = 5)
#> Error in attributes(.Data) <- c(attributes(.Data), attrib): length of 'dimnames' [2] not equal to array extent
```
