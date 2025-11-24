# Orientation tensor

3D orientation tensor, which characterize data distribution using
eigenvalue method. See (Watson 1966, Scheidegger 1965).

## Usage

``` r
is.ortensor(object)

as.ortensor(object)

ortensor(x, norm, w, shift)

# S3 method for class 'Vec3'
ortensor(x, norm = TRUE, w = NULL, shift = NULL)

# S3 method for class 'Line'
ortensor(x, norm = TRUE, w = NULL, shift = NULL)

# S3 method for class 'Ray'
ortensor(x, norm = TRUE, w = NULL, shift = NULL)

# S3 method for class 'Plane'
ortensor(x, norm = TRUE, w = NULL, shift = NULL)

# S3 method for class 'Pair'
ortensor(x, norm = TRUE, w = NULL, shift = FALSE)
```

## Arguments

- object:

  3x3 matrix

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- norm:

  logical. Whether the tensor should be normalized or not.

- w:

  numeric. weightings

- shift:

  logical. Only for `"Pair"` objects: Tensor is by default shifted
  towards positive eigenvalues, so it could be used as Scheidegger
  orientation tensor for plotting. When original Lisle tensor is needed,
  set shift to `FALSE`.

## Value

`ortensor()` returns an object of class `"ortensor"`

`is.ortensor` returns `TRUE` if `x` is an `"ortensor"` object, and
`FALSE` otherwise.

`as.ortensor` coerces a 3x3 matrix into an `"ortensor"` object.

## Details

The normalized orientation tensor is given as \$\$D = \frac{1}{n} (x_i,
y_i, z_i) (x_i, y_i, z_i)^T\$\$

## References

Watson, G. S. (1966). The Statistics of Orientation Data. The Journal of
Geology, 74(5), 786–797.

Scheidegger, A. E. (1964). The tectonic stress and tectonic motion
direction in Europe and Western Asia as calculated from earthquake fault
plane solutions. Bulletin of the Seismological Society of America,
54(5A), 1519–1528.
[doi:10.1785/BSSA05405A1519](https://doi.org/10.1785/BSSA05405A1519)

Lisle, R. (1989). The Statistical Analysis of Orthogonal Orientation
Data. The Journal of Geology, 97(3), 360-364.

## See also

[`inertia_tensor()`](https://tobiste.github.io/structr/reference/inertia.md)

Other ortensor:
[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md),
[`strain_shape`](https://tobiste.github.io/structr/reference/strain_shape.md)

## Examples

``` r
set.seed(20250411)
x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
ortensor(x, w = runif(nrow(x)))
#> Orientation tensor
#>            [,1]        [,2]        [,3]
#> [1,]  0.1953747 -0.22524425 -0.09751770
#> [2,] -0.2252442  0.69749852 -0.01237374
#> [3,] -0.0975177 -0.01237374  1.10031778

test <- as.ortensor(diag(3))
is.ortensor(test)
#> [1] TRUE

# Orientation tensor for Pairs
ortensor(angelier1990$TYM)
#> Orientation tensor
#>             [,1]        [,2]        [,3]
#> [1,] -0.33014781  0.17226487 -0.03483861
#> [2,]  0.17226487 -0.15036372 -0.06174581
#> [3,] -0.03483861 -0.06174581  0.48051152
```
