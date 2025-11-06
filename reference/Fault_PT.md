# Simple fault analysis

The PT-techniques is a graphical solution of the *Wallace-Bott
hypothesis*, i.e. fault slip occurs parallel to the maximum shear
stress. It calculates PT-axes, kinematic planes (also movement planes),
and the dihedra separation plane.

## Usage

``` r
Fault_PT(x, ptangle = 90)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- ptangle:

  numeric. angle between P and T axes in degrees (90° by default).

## Value

list. `p` and `t` are the P and T axes as `"Line"` objects, `m` and `d`
are the M-planes and the dihedra separation planes as `"Plane"` objects

## See also

[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)

## Examples

``` r
f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
Fault_PT(f)
#> $p
#> Line object (n = 3):
#>       azimuth   plunge
#> [1,] 314.9690 75.19665
#> [2,] 248.4545 15.32834
#> [3,] 342.4517 46.65113
#> 
#> $t
#> Line object (n = 3):
#>       azimuth   plunge
#> [1,] 116.2067 14.04868
#> [2,] 345.9919 25.60102
#> [3,] 241.3308 10.31893
#> 
#> $m
#> Plane object (n = 3):
#>      dip_direction      dip
#> [1,]      27.35344 85.42739
#> [2,]     310.64119 30.43222
#> [3,]     322.06010 48.49732
#> 
#> $d
#> Plane object (n = 3):
#>      dip_direction      dip
#> [1,]      289.7676 31.20884
#> [2,]      208.9071 83.18716
#> [3,]      210.2233 67.19866
#> 
```
