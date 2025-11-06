# Apparent dip direction

Apparent dip direction

## Usage

``` r
Plane_from_apparent_dips(a1, a2)
```

## Arguments

- a1, a2:

  `"Line"` objects containing the apparent dips and dip directions

## Value

`"Plane"` object

## Examples

``` r
a1 <- Line(45, 22)
a2 <- Line(352, 10)
res <- Plane_from_apparent_dips(a1, a2)

stereoplot()
points(rbind(a1, a2))
lines(res, lty = 2)
```
