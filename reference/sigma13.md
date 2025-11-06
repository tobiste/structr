# Principle stresses

calculates the magnitudes and directions of the principal stresses
\\\sigma_1\\ and \\\sigma_3\\

## Usage

``` r
sigma13(sigma_x, sigma_z, tau_xz)
```

## Arguments

- sigma_x:

  numeric. Magnitude of normal stress acting in the horizontal direction

- sigma_z:

  numeric. Magnitude of normal stress acting in the vertical direction

- tau_xz:

  numeric. Magnitude of shear stress acting on the same plane as
  `"sigma_x"`

## Value

A four-element list containing

- `"sigma1"`:

  magnitude of major principal stress

- `"sigma3"`:

  magnitude of minor principal stress

- `"theta1"`:

  direction of major principal stress (degrees)

- `"theta3"`:

  direction of minor principal stress (degrees)

## Author

Kyle Elmy and Jim Kaklamanos

## Examples

``` r
sigma13(sigma_x = 80, sigma_z = 120, tau_xz = 20)
#> $sigma1
#> [1] 128.2843
#> 
#> $sigma3
#> [1] 71.71573
#> 
#> $theta1
#> [1] 22.5
#> 
#> $theta3
#> [1] 112.5
#> 
```
