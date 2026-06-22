# Inertia tensor

Inertia tensor

## Usage

``` r
# S3 method for class 'spherical'
inertia_tensor(x, w = NULL)

inertia_tensor(x, w = NULL)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- w:

  numeric. weightings

## Value

3 x 3 matrix

## Details

\$\$D = n - (x_i, y_i, z_i) (x_i, y_i, z_i)^T\$\$

## See also

[`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md)

## Examples

``` r
set.seed(20250411)
x <- rfb(100, mu = Line(120, 50), k = 1, A = diag(c(10, 0, 0)))
inertia_tensor(x)
#>           [,1]      [,2]      [,3]
#> [1,]  91.20986 111.41893 102.75817
#> [2,] 111.41893  61.59032 104.72248
#> [3,] 102.75817 104.72248  47.19982
```
