# Structure classes

`Vec3`, `Line`, `Ray`, `Plane`, `"Pair"` and `Fault` create or convert a
`"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, and `"Fault"` S3 class
object, respectively, from the given set of values.

`as.Vec3`, `as.Line`, `as.Ray`, `as.Plane`, `as.Pair`, and `as.Fault`
attempt to coerce its argument into a `"Vec3"`, `"Line"`, `"Ray"`,
`"Plane"`, and `"Pair"`, and `"Fault"` S3 class object, respectively.

`is.Vec3`, `is.Line`, `is.Ray`, `is.Plane`, `is.Pair`, and `is.Fault`
test if its argument is a `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, and
`"Pair"`, and `"Fault"` S3 class object, respectively.

## Usage

``` r
is.spherical(x)

is.Vec3(x)

is.Line(x)

is.Ray(x)

is.Plane(x)

is.Pair(x)

is.Fault(x)

as.spherical(x, .class = NULL)

as.Vec3(x)

as.Line(x)

as.Ray(x)

as.Plane(x)

as.Pair(x)

as.Fault(x)

Vec3(x, y, z)

Line(x, plunge)

Ray(x, plunge, sense = NULL)

Plane(x, dip)

Pair(x, y, azimuth, plunge, correction = FALSE)

Fault(x, y, azimuth, plunge, sense, correction = FALSE)

Spherical(x, .class)
```

## Arguments

- x, y:

  object of class `"Line"`, `"Ray"`, `"Plane"`, and `"Pair"`, and
  `"Fault"` or numeric vector or array containing the spherical
  coordinates

- .class:

  character. Spherical class the object should be coerced to.

- sense:

  integer. Sense of the line (e.g.on a fault plane). Either `1`or `-1`
  for down (normal offset) or up (reverse offset), respectively. The
  "sense" is the sign of the fault's rake (see
  [`Fault_from_rake()`](https://tobiste.github.io/structr/reference/fault_from_rake.md)
  for details). Can also be a character with `"n"` (for normal) and
  `"r"` for "reverse".

- azimuth, plunge, z, dip:

  numeric vectors of the spherical coordinates

- correction:

  logical. If `TRUE` (default), both the fault plane and slip vector
  will be rotated so that the slip vector lies on the fault plane by
  minimizing the angle between the slip and the plane normal vector. See
  [`correct_pair()`](https://tobiste.github.io/structr/reference/pair_correct.md)
  for details.

## Details

`is.Vec3`, `is.Line`, `is.Plane`, `"is.Ray"`, `is.Pair`, and `is.Fault`
return `TRUE` if its arguments are an object of class `"Vec3"`,
`"Line"`, `"Ray"`, `"Plane"`, `"Pair"` or `"Fault"`, respectively, and
`FALSE` otherwise.

`is.spherical` returns `TRUE` if the argument's class is one of
`"Vec3()"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or `"Fault"` and
`FALSE` otherwise

`as.Vec3()`, `as.Line`, `as.Ray`, `as.Plane`, `as.Pair`, and `as.Fault`
are is generic functions.

A `Line` extends infinitely in both directions (equivalent to an axis in
2D), for example: principal stress directions, strain ellipsoid
directions (e.g. stretching lineation), intersection, fault striae,
crystallographic axes.

A `Ray` is a line with a single start point and extends indefinitely in
only one direction (equivalent to a direction in 2D): e.g. slip
direction, paleomagnetic direction (unless reversals are involved).

## Examples

``` r
x <- Line(120, 50) # create line
is.Line(x) # test if line
#> [1] TRUE
Plane(x) # convert to plane
#> Plane object (n = 1):
#> dip_direction           dip 
#>           300            40 
as.Plane(x) # assign as plane (note the difference to Pane(x))
#> Plane object (n = 1):
#> dip_direction           dip 
#>           120            50 

Pair(c(120, 120, 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23))
#> Pair object (n = 3):
#>      dip_direction dip azimuth plunge
#> [1,]           120  60     110     58
#> [2,]           120  60      25      9
#> [3,]           100  50      30     23
Fault(c("a" = 120, "b" = 120, "c" = 100), c(60, 60, 50), c(110, 25, 30), c(58, 9, 23), c(1, -1, 1))
#> Fault object (n = 3):
#>   dip_direction dip azimuth plunge sense
#> a           120  60     110     58     1
#> b           120  60      25      9    -1
#> c           100  50      30     23     1
```
