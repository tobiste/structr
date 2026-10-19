# Intersection of two planes

Intersection of two planes

## Usage

``` r
plane_intersection(x, y)
```

## Arguments

- x, y:

  Plane objects

## Value

Ray object representing the intersection of the two planes

## See also

[`crossprod.Plane()`](https://tobiste.github.io/structr/reference/vecmath.md)

## Examples

``` r
p1 <- Plane(120, 55)
p2 <- Plane(30, 45)
int <- plane_intersection(p1, p2)
print(int)
#> Ray object (n = 1):
#>   azimuth    plunge 
#> 245.00000 -39.32269 

stereoplot()
lines(rbind(p1, p2))
points(int, col = 'red', pch = 16)
```
