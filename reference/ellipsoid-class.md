# Ellipsoid class

In deformation analysis, the quadratic forms of the three-dimensional
stretches are represented by the `ellipsoid` class. It can be used to
represents either ellipsoid objects or finite strain ellipsoids.

## Usage

``` r
is.ellipsoid(x)

as.ellipsoid(x)

ellipsoid(x, left, ...)

# Default S3 method
ellipsoid(x, left = NULL, ...)

# S3 method for class 'defgrad'
ellipsoid(x, left = TRUE, ...)
```

## Arguments

- x:

  Either a matrix or a `"defgrad"` object

- left:

  logical. `TRUE` for left Cauchy–Green deformation tensor (or Finger
  deformation tensor, the default), `FALSE` for for right Cauchy–Green
  deformation tensor (or Green’s deformation tensor)

- ...:

  optional parameters passed to function call.

## Value

`is.ellipsoid` returns `TRUE` if `x` is an `"ellipsoid"` object, and
`FALSE` otherwise.

`as.ellipsoid` coerces a 3x3 matrix into an `"ellipsoid"` object.

## Details

The eigenvalues \\\lambda\\ of the deformation matrix are the quadratic
forms of the principal stretches \\s\\ (\\s = 1 + \epsilon = l/l_0\\).

## See also

[`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md)

Other ellipsoid:
[`ellipsoid-params`](https://tobiste.github.io/structr/reference/ellipsoid-params.md),
[`ellipsoid_from_stretch()`](https://tobiste.github.io/structr/reference/ellipsoid_from_stretch.md)

## Examples

``` r
test <- as.ellipsoid(diag(3))
is.ellipsoid(test)
#> [1] TRUE

R <- defgrad_from_ratio(2, 3)
ellipsoid(R)
#> Ellipsoid tensor
#>          [,1]     [,2]      [,3]
#> [1,] 5.241483 0.000000 0.0000000
#> [2,] 0.000000 1.310371 0.0000000
#> [3,] 0.000000 0.000000 0.1455967
```
