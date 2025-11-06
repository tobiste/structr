# Extract components of fault object

`Fault_plane()` extracts the orientation of the fault plane,
`Fault_slip()` extracts the orientation of the slip vector, and
`Fault_rake()` extracts the rake of the fault, i.e. the angle between
fault slip vector and fault strike (measured clockwise from strike!)
`Fault_sense()` extracts the fault sense from the rake (1: normal, -1:
reverse)

## Usage

``` r
Fault_rake(x)

Pair_plane(x)

Fault_slip(x)

Fault_plane(x)

Fault_sense(x, steps = 8)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- steps:

  Integer. Either 2, 4, or 8 steps for parsing the fault sense

## Value

numeric. `"Plane"`, `"Ray"` or angle in degrees, respectively

## Examples

``` r
f <- Fault(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
Fault_plane(f)
#> Plane object (n = 3):
#>      dip_direction dip
#> [1,]           120  60
#> [2,]           120  60
#> [3,]           100  50
Fault_slip(f)
#> Line object (n = 3):
#>      azimuth plunge
#> [1,]     110     58
#> [2,]      25      9
#> [3,]      30     23
Fault_rake(f)
#> [1]  84.72020 -10.28562  30.11825
Fault_sense(f, 2)
#> [1] Normal  Reverse Normal 
#> Levels: Normal Reverse
Fault_sense(f, 4)
#> [1] Normal-sinistral  Reverse-sinistral Normal-sinistral 
#> 4 Levels: Normal-sinistral Normal-dextral ... Reverse-sinistral
Fault_sense(f, 8)
#> [1] Normal-sinistral  Reverse-sinistral Sinistral        
#> 8 Levels: Sinistral Normal-sinistral Normal Normal-dextral ... Reverse-sinistral
```
