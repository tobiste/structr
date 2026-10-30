# Axial Cross Moments

Axial Cross Moments

## Usage

``` r
axial_cross_moments(x, nu)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- nu:

  Three Eigenvectors \\\nu\\ (object of class `Spherical`). See
  [`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)
  for details

## Value

3 x 3 matrix

## Details

Matrix of unit vectors `x` for Eigenvectors `\nu` is \$\$m\_{ij} =
\frac{1}{n} \sum\_{k=1}^{n} \left( \mathbf{x}\_k \cdot \mathbf{\nu}\_i
\right)^2 \left( \mathbf{x}\_k \cdot \mathbf{\nu}\_j \right)^2 , \qquad
i, j = 1, 2, 3 \$\$

## Examples

``` r
x <- example_planes
axial_cross_moments(x, ot_eigen(x)$vectors)
#>            [,1]       [,2]       [,3]
#> [1,] 0.48997614 0.09593684 0.04384632
#> [2,] 0.09593684 0.17415304 0.01803057
#> [3,] 0.04384632 0.01803057 0.02024336
```
