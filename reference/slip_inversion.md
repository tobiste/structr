# Stress Inversion for Fault-Slip Data

Convenient function for linear and direction inversion of fault-slip
data to derive the reduced stress tensor.

## Usage

``` r
slip_inversion(x, method = c("michael", "angelier"), ...)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- method:

  character. The inversion algorithm. One of `"michael"` (the default)
  for a bootstrapped linear inversion, or `"angelier"` for a iterative
  direct inversion.

- ...:

  arguments passed to
  [`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md)
  or
  [`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md)

## Value

list. See
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md)
and
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md)

## References

Angelier, J. (1990). Inversion of field data in fault tectonics to
obtain the regional stress—III. A new rapid direct inversion method by
analytical means. Geophys. J. Int, 103, 363–376.
<https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

Michael, A. J. (1984). Determination of stress from slip data: Faults
and folds. Journal of Geophysical Research: Solid Earth, 89(B13),
11517–11526. <https://doi.org/10.1029/JB089iB13p11517>

## See also

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md)

## Examples

``` r
set.seed(20250411)
# Use Angelier examples
par(mfrow = c(1, length(angelier1990)))
invisible(lapply(angelier1990, function(x){

# Linear inversion
TYM_michael <- slip_inversion(x, method = "michael")

# Direct inversion
TYM_angelier <- slip_inversion(x, method = "angelier")

stereoplot(guides = FALSE)
fault_plot(x, col = "gray80")
points(TYM_michael$principal_axes, pch = 3, col = 1:3)
points(TYM_angelier$principal_axes, pch = 2, col = 1:3)
legend("topleft", col = 1, legend = c("Michael (1984)", "Angelier (1990)"), pch = c(3, 2))
}))
```
