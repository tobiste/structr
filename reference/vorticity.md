# Flow Apophyses, Vorticity, and Instantaneous Stretching Axes

Computes flow apophyses, vorticity vectors, kinematic vorticity numbers,
and instantaneous stretching axes from eigenvalues and eigenvectors of
velocity gradient tensor

## Usage

``` r
flow_apophyses(x)

vorticity_axis(x)

kinematic_vorticity_from_velgrad(x)

instantaneous_stretching(x)

instantaneous_stretching_axes(x)
```

## Arguments

- x:

  object of class `"velgrad"`

## Value

`vorticity_axis` and `instantaneous_stretching_axes` return `"Vec3"`
object; `"flow_apophyses"` returns `"Plane"` object;
`kinematic_vorticity_from_velgrad` and `instantaneous_stretching` return
numeric.

## Examples

``` r
d <- defgrad_from_generalshear(k = 2.5, gamma = 0.9)
l <- velgrad(d, time = 10)

flow_apophyses(l)
#> Plane object (n = 2):
#>   dip_direction dip
#> y      135.5139  90
#> y      270.0000  90
vorticity_axis(l)
#> Vector (Vec3) object (n = 1):
#> x y z 
#> 0 0 1 
kinematic_vorticity_from_velgrad(l)
#> [1] 0.7007364
instantaneous_stretching_axes(l)
#> Vector (Vec3) object (n = 3):
#>              x          y z
#> [1,] 0.3785365  0.9255864 0
#> [2,] 0.9255864 -0.3785365 0
#> [3,] 0.0000000  0.0000000 1
instantaneous_stretching(l)
#> [1]  0.11003269 -0.01840362 -0.09162907
```
