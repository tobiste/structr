# Vector Rotation

Vector Rotation

## Usage

``` r
rotate(x, rotaxis, rotangle)

# S3 method for class 'Vec3'
rotate(x, rotaxis, rotangle)

# S3 method for class 'Ray'
rotate(x, rotaxis, rotangle)

# S3 method for class 'Line'
rotate(x, rotaxis, rotangle)

# S3 method for class 'Plane'
rotate(x, rotaxis, rotangle)

# S3 method for class 'Pair'
rotate(x, rotaxis, rotangle)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- rotaxis:

  Axis of rotation given as object of class `"Vec3"`, `"Line"`, `"Ray"`,
  or `"Plane"`.

- rotangle:

  Angle of rotation in radians for `"Vec3"` objects and in degrees for
  `"Line"`, `"Ray"` and `"Plane"` objects.

## Value

objects of same class as `x`

## See also

[`defgrad_from_axisangle()`](https://tobiste.github.io/structr/reference/defgrad.md);
[`transform_linear()`](https://tobiste.github.io/structr/reference/vecmath.md)

## Examples

``` r
vec1 <- Vec3(1, 0, 0)
vec2 <- Vec3(0, 0, 1)
rotate(vec1, vec2, pi / 2)
#> Vector (Vec3) object (n = 1):
#>            x            y            z 
#> 6.123234e-17 1.000000e+00 0.000000e+00 

# rotate Fault data (sense of motion changes!)
rotate(simongomez[1:5, ], Ray(90, 10), 80)
#> Fault object (n = 5):
#>      dip_direction      dip  azimuth   plunge sense
#> [1,]     250.46711 42.09286 250.9936 42.09166    -1
#> [2,]      92.59349 82.28746 142.7678 78.06032    -1
#> [3,]     281.39043 44.48230 249.1242 39.70735    -1
#> [4,]     293.63164 59.35347 246.4979 48.94592     1
#> [5,]     296.64621 46.49228 259.0021 39.83414     1
```
