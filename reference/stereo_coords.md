# Stereographic and Equal-Area Projection

Transformation of spherical coordinates into the stereographic or
equal-area projection

## Usage

``` r
stereo_coords(az, inc, upper.hem = NULL, earea = NULL, radius = 1)
```

## Arguments

- az, inc:

  numeric vectors. Azimuth and Inclination in degrees.

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`). Defaults to
  `getOption("structr.upper.hem")`.

- earea:

  logical. Projection, either `TRUE` for Lambert equal-area projection,
  or `FALSE` for meridional stereographic projection. Defaults to
  `getOption("structr.earea")`.

- radius:

  numeric. Radius of circle. Defaults to `getOption("structr.radius")`.

## Value

two-column vector with the transformed coordinates

## See also

[structr-options](https://tobiste.github.io/structr/reference/structr-options.md)

## Examples

``` r
stereo_coords(90, 10)
#>            x            y
#> inc 0.909039 5.566258e-17
stereo_coords(90, 10, earea = TRUE, upper.hem = TRUE)
#>             x             y
#> inc -0.909039 -1.669877e-16
```
