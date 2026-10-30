# Angle between vectors

Angle between vectors

## Usage

``` r
angle(x, y, axial, ...)

# S3 method for class 'spherical'
angle(x, y, axial = NULL, ...)
```

## Arguments

- x, y:

  objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- axial:

  logical. Whether the angle is measured between axial (`Line`, `Plane`)
  or vectorial (`Ray`, `Vec3`) data. Axial angles are folded into the
  range 0–90°, vectorial angles range from 0–180°. For spherical
  objects, the default (`NULL`) determines this from the classes of `x`
  and `y`: the angle is axial if either is a `Line` or `Plane`.

- ...:

  not used.

## Value

numeric. Angle in degrees if either `x` or `y` is of class `Line`,
`Ray`, or `Plane`, otherwise in radians. For a `Plane`, the angle is
measured to its pole.

## Examples

``` r
angle(Ray(0, 0), Ray(180, 0))
#> [1] 180
angle(Line(0, 0), Line(180, 0))
#> [1] 0

# force vectorial treatment of lines:
angle(Line(0, 0), Line(180, 0), axial = FALSE)
#> [1] 180
```
