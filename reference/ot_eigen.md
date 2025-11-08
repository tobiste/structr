# Eigenvalues and Eigenvectors of a Set of Vectors

Decomposition of Orientation Tensor Eigenvectors and Eigenvalues

## Usage

``` r
ot_eigen(x, scaled = FALSE, ...)

# S3 method for class 'spherical'
ot_eigen(x, scaled = FALSE, ...)

# S3 method for class 'ortensor'
ot_eigen(x, scaled = FALSE, ...)
```

## Arguments

- x:

  either an object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`,
  `"Pair"`, or `"Fault"` where the rows are the observations and the
  columns are the coordinates, or an `"ortensor"` object.

- scaled:

  logical. Whether the Eigenvectors should be scaled by the Eigenvalues
  (only effective if `x` is in Cartesian coordinates).

- ...:

  additional arguments passed to
  [`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md)
  (ignored if `x` is `"ortensor"` object).

## Value

list containing

- `values`:

  Eigenvalues

- `vectors`:

  Eigenvectors in coordinate system of `x`

## See also

Other ortensor:
[`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md),
[`strain_shape`](https://tobiste.github.io/structr/reference/strain_shape.md)

## Examples

``` r
set.seed(20250411)
mu <- rvmf(n = 1)
x <- rfb(100, mu = mu, k = 1, A = diag(c(10, 0, 0)))
x_eigen <- ot_eigen(x)
x_eigen
#> eigen() decomposition
#> $values
#> [1] 0.48995265 0.46558493 0.04446242
#> 
#> $vectors
#> Vector (Vec3) object (n = 3):
#>               x          y         z
#> [1,] -0.6917239  0.1673371 0.7025072
#> [2,]  0.6572339  0.5490008 0.5163737
#> [3,]  0.2992685 -0.8188996 0.4897365
#> 
plot(x, col = "grey")
points(mu, col = 4)
text(mu, labels = "Mean", col = 4, pos = 4)
points(x_eigen$vectors, col = c(1, 2, 3))
text(x_eigen$vectors, col = c(1, 2, 3), labels = c("E1", "E2", "E3"), pos = 4)
```
