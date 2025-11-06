# Principal stresses from 2D stress components

Determines the principal stresses and their orientations from the stress
components \$\$sigma_x\$\$, \$\$\sigma_y\$\$, \$\$\tau\_{xy}\$\$.

## Usage

``` r
PR_stress(sigma_x, sigma_y, tau_xy)
```

## Arguments

- sigma_x:

  numeric. Magnitude of normal stress acting in the horizontal direction

- sigma_y:

  numeric. Magnitude of normal stress acting on plane facing in Y
  direction (\$\$\sigma_y\$\$).

- tau_xy:

  numeric. Magnitude of shear stress acting on planes facing X and Y
  (\$\$\tau\_{xy}\$\$).

## Value

angle in degrees

## References

Richard J. Lisle (1999)

## Examples

``` r
PR_stress(sigma_x = 80, sigma_y = 120, tau_xy = 20)
#>     sigma1   sigma2 theta
#> 1 128.2843 71.71573  67.5
```
