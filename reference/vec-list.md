# List of vectors

Creates a list of Cartesian vectors from an spherical objects. This is a
convenience function to link with the package `geologyGeometry` by J. R.
Davis

## Usage

``` r
vec_list(x)

list_vec(ls)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- ls:

  list of 3-element vectors (Cartesian coordinates x, y, z)

## Value

`vec_list` returns a list of 3-element vectors (Cartesian coordinates x,
y, z). `list_vec` returns a `"Vec3"` object

## Examples

``` r
ls <- vec_list(example_lines[1:5, ])
print(ls)
#> [[1]]
#>         x         y         z 
#> 0.5727204 0.7882819 0.2249511 
#> 
#> [[2]]
#>         x         y         z 
#> 0.4682901 0.8448178 0.2588190 
#> 
#> [[3]]
#>         x         y         z 
#> 0.2674497 0.9327081 0.2419219 
#> 
#> [[4]]
#>         x         y         z 
#> 0.1641876 0.9311540 0.3255682 
#> 
#> [[5]]
#>         x         y         z 
#> 0.4341533 0.8520738 0.2923717 
#> 

list_vec(ls)
#> Vector (Vec3) object (n = 5):
#>              x         y         z
#> [1,] 0.5727204 0.7882819 0.2249511
#> [2,] 0.4682901 0.8448178 0.2588190
#> [3,] 0.2674497 0.9327081 0.2419219
#> [4,] 0.1641876 0.9311540 0.3255682
#> [5,] 0.4341533 0.8520738 0.2923717
```
