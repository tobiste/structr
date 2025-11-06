# Fault from plane and rake

Fault from plane and rake

## Usage

``` r
Fault_from_rake(p, rake, sense = NULL, ...)
```

## Arguments

- p:

  object of class `"Plane"`

- rake:

  Angle (in degrees) between fault strike and lineation. Measured
  clockwise from the strike, i.e. down is positive; values between 0 and
  360° (or −180° and 180°)

- sense:

  Either 1 (for normal fault movement) or -1 (reverse fault movement).
  Use this only if `rake` values are in range of 0 and 180°.

- ...:

  optional arguments passed to
  [`Fault()`](https://tobiste.github.io/structr/reference/classes.md)

## Value

`"Fault"` object

## Details

Rake is used to describe the direction of fault motion with respect to
the strike.

Measured clockwise from the strike, down is positive; values between 0
and 360° (or −180° and 180°):

- left-lateral strike slip: rake near 0°

- right-lateral strike slip: rake near 180°

- normal: rake near 90°

- reverse/thrust: rake near -90° (270°)

## Examples

``` r
fr <- Fault_from_rake(Plane(c(120, 120, 100, 0), c(60, 60, 50, 40)),
  rake = c(84.7202, -10, 30, 180)
)
plot(fr, col = 1:4)
legend("topleft", legend = Fault_sense(fr, 8), col = 1:4, pch = 16)


fr2 <- Fault_from_rake(Plane(c(90, 90, 90), c(80, 40, 10)),
  rake = c(10, 20, 90), sense = c(1, 1, -1)
)
plot(fr2, col = 1:3)
legend("topleft", legend = Fault_sense(fr2, 8), col = 1:3, pch = 16)
```
