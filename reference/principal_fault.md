# Calculate Principal Fault planes from stress vectors and friction

Calculate Principal Fault planes from stress vectors and friction

## Usage

``` r
principal_fault(s1, s3, friction = 0.6)
```

## Arguments

- s1, s3:

  Principal stress vectors as spherical objects

- friction:

  numeric. Coefficient of friction (0.6 by default)

## Value

`"Fault"` object

## See also

[`stereo_fms()`](https://tobiste.github.io/structr/reference/stereo_fms.md)
to plot the focal mechanism solution

## Examples

``` r
s <- reduced_stress(angelier1990$TYM)
stress <- sigma2stress(s)
pr_TYM <- principal_fault(stress$axes[1, ], stress$axes[3, ])

stereoplot()
fault_plot(angelier1990$TYM, col = "grey")
fault_plot(pr_TYM, col = "red")
points(stress$axes, pch = 16)


# Focal mechanism plot
stereoplot()
stereo_fms(s)
fault_plot(angelier1990$TYM, col = "grey30", type = "hoeppener")
points(stress$axes, col = 2:4, pch = 16)
```
