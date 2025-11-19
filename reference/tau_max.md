# Maximum in-plane shear stress

calculates the magnitude and direction of the maximum in-plane shear
stress

## Usage

``` r
tau_max(sigma_x, sigma_z, tau_xz)
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

A two-element list containing

- `"tauMax"`:

  maximum in-plane shear stress

- `"theta"`:

  angle of maximum in-plane shear stress (in degrees)

## Author

Kyle Elmy and Jim Kaklamanos

## Examples

``` r
tau_max(sigma_x = 80, sigma_z = 120, tau_xz = 20)
#> $tau_max
#> function (sigma_x, sigma_z, tau_xz) 
#> {
#>     tauMax <- sqrt((((sigma_z - sigma_x)/2)^2) + ((tau_xz)^2))
#>     x <- (1/(1 + ((2 * tau_xz)/(sigma_z - sigma_x))^2))
#>     theta1 <- acosd(sqrt(x))/2
#>     theta <- theta1 + 45
#>     return(list(tau_max = tau_max, theta = theta))
#> }
#> <bytecode: 0x562595b174b0>
#> <environment: namespace:structr>
#> 
#> $theta
#> [1] 67.5
#> 
```
