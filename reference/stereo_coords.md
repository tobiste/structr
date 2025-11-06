# Stereographic Projection

Transformation of spherical coordinates into the stereographic
projection

## Usage

``` r
stereo_coords(az, inc, upper.hem = FALSE, earea = TRUE, r = 1)
```

## Arguments

- az, inc:

  numeric vectors. Azimuth and Inclination in degrees.

- upper.hem:

  logical. Whether the projection is shown for upper hemisphere (`TRUE`)
  or lower hemisphere (`FALSE`, the default).

- earea:

  logical `TRUE` for Lambert equal-area projection (also "Schmidt net";
  the default), or `FALSE` for meridional stereographic projection (also
  "Wulff net" or "Stereonet").

- r:

  numeric. Radius of circle. Default is `1` for unit circle.

## Value

two-column vector with the transformed coordinates

## Examples

``` r
stereo_coords(90, 10)
#>            x            y
#> inc 0.909039 5.566258e-17
stereo_coords(90, 10, earea = TRUE, upper.hem = TRUE)
#>             x             y
#> inc -0.909039 -1.669877e-16
```
