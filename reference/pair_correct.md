# Orthogonalization of plane and line measurement

Both the line and the plane are rotated in opposite directions by half
the angle between the slip and the plane normal vector about the vector
that is perpendicular to both.

## Usage

``` r
misfit_pair(x)

correct_pair(x)
```

## Arguments

- x:

  object of class `"Pair"` or `"Fault"`

## Value

`misfit_pair` returns a list with the orthogonalized plane and line
measurements (as 3d vectors) and the misfit angles (in radians).
`correct_pair` returns the orthogonalized vectors.

## Examples

``` r
p <- Pair(120, 60, 110, 58, correction = FALSE)
misfit_pair(p)
#> $fvec
#> Vector (Vec3) object (n = 1):
#>          x          y          z 
#>  0.4306083 -0.7432653  0.5119893 
#> 
#> $lvec
#> Vector (Vec3) object (n = 1):
#>          x          y          z 
#> -0.1752490  0.4876331  0.8552787 
#> 
#> $misfit
#> [1] 0.02793105
#> 

correct_pair(p)
#> Pair object (n = 1):
#> dip_direction           dip       azimuth        plunge 
#>     120.08572      59.20357     109.76778      58.79054 
```
