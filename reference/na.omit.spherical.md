# Handle Missing Values in Spherical Objects

Handle Missing Values in Spherical Objects

## Usage

``` r
# S3 method for class 'spherical'
na.omit(object, ...)
```

## Arguments

- object:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`.

- ...:

  further arguments special methods could require.

## Examples

``` r
x <- Line(c(120, NA, 100), c(50, 60, 70))
na.omit(x)
#> Line object (n = 2):
#>      azimuth plunge
#> [1,]     120     50
#> [2,]     100     70
```
