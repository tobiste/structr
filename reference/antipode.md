# Antipode vector

Returns the opposite vector

## Usage

``` r
antipode(x, ...)

# Default S3 method
antipode(x, ...)

# S3 method for class 'Vec3'
antipode(x, ...)

# S3 method for class 'Line'
antipode(x, ...)

# S3 method for class 'Ray'
antipode(x, ...)

# S3 method for class 'Plane'
antipode(x, ...)

# S3 method for class 'Pair'
antipode(x, ...)

# S3 method for class 'Fault'
antipode(x, ...)
```

## Arguments

- x:

  a spherical object

- ...:

  arguments passed to function call

## Examples

``` r
Line(120, 55) |> antipode()
#> Line object (n = 1):
#> azimuth  plunge 
#>     120      55 
Ray(120, 55) |> antipode()
#> Ray object (n = 1):
#> azimuth  plunge 
#>     300     -55 
```
