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
#>              x          y           z
#> [1,] 0.9428027  0.2728321  0.19153539
#> [2,] 0.6241873  0.7194017 -0.30471512
#> [3,] 0.9562071  0.2774389  0.09325005
#> [4,] 0.9110783 -0.2832872  0.29947405
#> [5,] 0.7246633  0.3258554 -0.60719138
#> [6,] 0.9216500 -0.2600964 -0.28794290
tail(x)
#> Vector (Vec3) object (n = 6):
#>              x          y            z
#> [1,] 0.7246633  0.3258554 -0.607191382
#> [2,] 0.9216500 -0.2600964 -0.287942903
#> [3,] 0.4542605  0.4783322 -0.751562168
#> [4,] 0.9680807  0.2505073 -0.008117167
#> [5,] 0.9400175 -0.3357350  0.060408045
#> [6,] 0.6450888  0.2964477 -0.704257908
```
