# Stress transformation

calculates the magnitudes of the normal stress and the shear stress

## Usage

``` r
stress_transformation(
  theta,
  sigma_x = NA,
  sigma_z = NA,
  tau_xz = NA,
  sigma1 = NA,
  sigma3 = NA
)
```

## Arguments

- theta:

  numeric. Angles (degrees); defaults to 0-180 in increments of 1

- sigma_x:

  numeric. Magnitude of normal stress acting in the horizontal direction

- sigma_z:

  numeric. Magnitude of normal stress acting in the vertical direction

- tau_xz:

  numeric. Magnitude of shear stress acting on the same plane as
  `"sigma_x"`

- sigma1:

  numeric. Magnitude of major principal stress

- sigma3:

  numeric. Magnitude of minor principal stress

## Value

A two-element list containing

- `"normal"`:

  normal stress on an inclined plane

- `"shear"`:

  shear stress on an inclined plane

## Note

In addition to theta, One of the following two sets of data must be
entered:

1.  `"sigma_x"`, `"sigma_z"`, `"tau_xz"`

2.  `"sigma1"`, `"sigma3"`

If theta is entered in conjunction with `"sigma_x"`, `"sigma_z"`, and
`"tau_xz"`, it is interpreted as the angle of inclination above the
horizontal. If theta is entered in conjunction with the principal
stresses, then it is interpreted as the angle of inclination above the
major principal plane.

## Author

Kyle Elmy and Jim Kaklamanos

## Examples

``` r
stress_transformation(sigma_x = 80, sigma_z = 120, tau_xz = 20, theta = 78)
#> $normal
#> [1] 89.86382
#> 
#> $shear
#> [1] 26.40564
#> 
```
