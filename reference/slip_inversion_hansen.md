# 9D Direct Inversion of Fault Slip

Direct inversion of stress, strain or strain rate including vorticity
using 9D parameter space using the method by Hansen (2013)

## Usage

``` r
slip_inversion_hansen(x, flip = FALSE)
```

## Arguments

- x:

  object of class `"Pair"` or `"Fault"`

- flip:

  logical. Flip if you want to have the negative stress tensor, i.e.
  sigma 1 and 3 will be flipped.

## Value

list

## References

Hansen, J. A. (2013). Direct inversion of stress, strain or strain rate
including vorticity: A linear method of homogenous fault-slip data
inversion independent of adopted hypothesis. Journal of Structural
Geology, 51, 3–13. https://doi.org/10.1016/j.jsg.2013.03.014

## See also

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md)

## Examples

``` r
par(mfrow = c(1, 2))
invisible(lapply(angelier1990, function(x){

res <- slip_inversion_hansen(x, FALSE)

plot(x, col = 'lightgrey')
title(sub = paste0("R = ", round(res$phi, 2)))
points(res$principal_axes, col = 2:4, pch = 16, cex = 2)
text(res$principal_axes, labels = 1:3, col = 2:4, adj = -1)
points(res$vorticity_axis, col = 5, pch = 17, cex = 2)
}))
#> Error in slip_inversion_hansen(x, FALSE): object 'MSort' not found
```
