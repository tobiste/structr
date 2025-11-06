# Return the First or Last Parts of an Object

Returns the first or last parts of a vector.

## Usage

``` r
# S3 method for class 'spherical'
head(x, n = 6L, ...)

# S3 method for class 'spherical'
tail(x, n = 6L, ...)
```

## Arguments

- x:

  objects of class `"Vec3"`, `"Line"`, `"Plane"`, `"Pair"`, or `"Fault`

- n:

  an integer vector of length up to `dim(x)` (or 1, for non-dimensioned
  objects). A `logical` is silently coerced to integer. Values specify
  the indices to be selected in the corresponding dimension (or along
  the length) of the object. A positive value of `n[i]` includes the
  first/last `n[i]` indices in that dimension, while a negative value
  excludes the last/first `abs(n[i])`, including all remaining indices.
  `NA` or non-specified values (when `length(n) < length(dim(x))`)
  select all indices in that dimension. Must contain at least one
  non-missing value.

- ...:

  arguments to be passed to or from other methods.

## Examples

``` r
x <- rvmf(n = 10)
head(x)
#> Vector (Vec3) object (n = 6):
#>              x          y           z
#> [1,] 0.7303984 0.28941147 -0.61867539
#> [2,] 0.8814433 0.39564139 -0.25792560
#> [3,] 0.9146064 0.10018765  0.39173650
#> [4,] 0.5910464 0.09677367  0.80081144
#> [5,] 0.8830778 0.45301240  0.12228394
#> [6,] 0.8021508 0.59634824  0.03037924
tail(x)
#> Vector (Vec3) object (n = 6):
#>              x           y           z
#> [1,] 0.8830778  0.45301240  0.12228394
#> [2,] 0.8021508  0.59634824  0.03037924
#> [3,] 0.8474214  0.03356496  0.52985887
#> [4,] 0.9276049  0.31544268 -0.20011264
#> [5,] 0.7842683 -0.44050598  0.43689552
#> [6,] 0.8932390  0.20969618 -0.39768281
```
