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
#> [1,]  89.99684 111.33511 104.70367
#> [2,] 111.33511  61.14329 102.89073
#> [3,] 104.70367 102.89073  48.85986
```
