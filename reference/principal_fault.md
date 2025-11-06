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

  numeric. Coefficient of friction

## Value

`"Fault"` object

## Examples

``` r
res_TYM <- slip_inversion(angelier1990$TYM, boot = 10)
pr_TYM <- principal_fault(res_TYM$principal_axes[1, ], res_TYM$principal_axes[3, ])

stereoplot()
fault_plot(angelier1990$TYM, col = "grey")
fault_plot(pr_TYM, col = "red")
points(res_TYM$principal_axes, pch = 16)
```
