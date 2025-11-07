# Fault from rake and quadrant notation

Fault from rake and quadrant notation

## Usage

``` r
Fault_from_rake_quadrant(
  p,
  rake,
  quadrant,
  type = c("plunge", "rake"),
  sense = NULL
)
```

## Arguments

- p:

  `"Plane"` object

- rake:

  numeric. Rake angle in degrees

- quadrant:

  character. Quadrant of plunge or rake direction

- type:

  character. Either `"plunge"` or `"rake"` for specifying which quadrant
  convention is used

- sense:

  Either 1 (for normal fault movement) or -1 (reverse fault movement).
  Only used when `type=="rake"`

## Value

`"Fault"` object

## Details

`type=="plunge"` This is the angle measured in the fault plane between
the strike given by either right or left-hand rule and the lineation.
The angle is recorded in a clockwise sense (looking down upon the fault
plane) and has a range from 0 to 180%deg;. The quadrant of plunge
indicates the direction of the strike from which the rake angle is
measured. The `sense` argument is ignored, as it it implied by the sign
of the rake.

`type=="rake"` Rake is the **acute** angle measured in the fault plane
between the strike of the fault and the lineation. Starting from the
strike line, the angle is measured in a sense which is down the dip of
the plane. Quadrant of rake indicate the direction of the strike from
which the rake angle is measured, i.e. whether right-hand or
left-handrule is followed. Angle ranges from 0 to 90 °. Use `sense`
argument to specify the sense of motion.

## See also

Other parse-orientations:
[`azimuth_to_cardinal()`](https://tobiste.github.io/structr/reference/azimuth_to_cardinal.md),
[`quadrant2dd()`](https://tobiste.github.io/structr/reference/quadrant2dd.md),
[`sense_from_strikeslip()`](https://tobiste.github.io/structr/reference/sense_from_strikeslip.md),
[`split()`](https://tobiste.github.io/structr/reference/split.md),
[`strikeslip_kinematics()`](https://tobiste.github.io/structr/reference/strikeslip_kinematics.md)

## Examples

``` r
dip <- c(5, 10, 15, 30, 40, 55, 65, 75, 90)
dip_dir <- c(180, 225, 270, 315, 360, 0, 45, 90, 135)
rake1 <- c(0, 45, 90, 135, 180, 45, 90, 135, 180)
plunge_quadrant <- c("E", "S", "W", "N", "E", "W", "E", "S", "W")
Fault_from_rake_quadrant(Plane(dip_dir, dip), rake1, plunge_quadrant, type = "plunge")
#> Fault object (n = 9):
#>       dip_direction dip    azimuth       plunge sense
#>  [1,]           180   5  90.000000 0.000000e+00     1
#>  [2,]           225  10 179.561451 7.053022e+00     1
#>  [3,]           270  15 270.000000 1.500000e+01     1
#>  [4,]           315  30   4.106605 2.070481e+01     1
#>  [5,]             0  40  90.000000 1.487542e-14     1
#>  [6,]             0  55 299.837566 3.539626e+01     1
#>  [7,]            45  65  45.000000 6.500000e+01     1
#>  [8,]            90  75 165.489181 4.307952e+01     1
#>  [9,]           135  90 225.000000 7.016709e-15     1

rake2 <- c(0, 45, 90, 45, 0, 45, 90, 45, 0)
rake_quadrant <- c("E", "S", "S", "E", "E", "W", "N", "S", "W")
Fault_from_rake_quadrant(Plane(dip_dir, dip), rake2, rake_quadrant, type = "rake")
#> Fault object (n = 9):
#>       dip_direction dip    azimuth    plunge sense
#>  [1,]           180   5  90.000000  0.000000     1
#>  [2,]           225  10 179.561451  7.053022     1
#>  [3,]           270  15 270.000000 15.000000     1
#>  [4,]           315  30   4.106605 20.704811     1
#>  [5,]             0  40 270.000000  0.000000     1
#>  [6,]             0  55 299.837566 35.396260     1
#>  [7,]            45  65  45.000000 65.000000     1
#>  [8,]            90  75 165.489181 43.079517     1
#>  [9,]           135  90  45.000000  0.000000     1
```
