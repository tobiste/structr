# Return the First or Last Parts of a Spherical Object

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
#>              x           y           z
#> [1,] 0.9194460  0.07680528 -0.38564221
#> [2,] 0.8944839  0.29503629 -0.33593484
#> [3,] 0.7650453  0.63719963  0.09317877
#> [4,] 0.6598253  0.06208816  0.74884950
#> [5,] 0.9500657 -0.07762652 -0.30224032
#> [6,] 0.8170436 -0.03978027  0.57520197
tail(x)
#> Vector (Vec3) object (n = 6):
#>              x           y           z
#> [1,] 0.9500657 -0.07762652 -0.30224032
#> [2,] 0.8170436 -0.03978027  0.57520197
#> [3,] 0.6663839  0.74429159  0.04430004
#> [4,] 0.8913334  0.15888109 -0.42459586
#> [5,] 0.8677729 -0.25526040 -0.42639448
#> [6,] 0.9878350 -0.11979314  0.09915423
```
