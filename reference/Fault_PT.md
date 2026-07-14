# Simple fault analysis

The PT-techniques is a graphical solution of the *Wallace-Bott
hypothesis*, i.e. fault slip occurs parallel to the maximum shear
stress. It calculates PT-axes, kinematic planes (also movement planes),
and the dihedra separation plane.

## Usage

``` r
Fault_PT(x, ptangle = 90)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates. Object must be complete, i.e. no `NA` values. For
  Michael's, Angelier's, and Yamaji-Sato's methods, at least 4 rows of
  fault measurements are required, while Hansen's method requires at
  least 7.

- ptangle:

  numeric. angle between P and T axes in degrees (90° by default).

## Value

list. `p` and `t` are the P and T axes as `"Line"` objects, `m` and `d`
are the M-planes and the dihedra separation planes as `"Plane"` objects

## See also

Other stress-inversion:
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_hansen_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen_boot.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md),
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md),
[`slip_inversion_yamaji_sato()`](https://tobiste.github.io/structr/reference/slip_inversion_yamaji_sato.md)

## Examples

``` r
par(mfrow = c(1, length(angelier1990)))
invisible(lapply(angelier1990, function(x) {
  xpt <- Fault_PT(x)

  stereoplot(guides = FALSE)
  angelier(x, col = "grey")
  points(xpt$p, pch = 16, cex = 0.6, col = 1)
  points(xpt$t, pch = 16, cex = 0.6, col = 2)
  stereo_confidence(xpt$p, pch = 16, cex = 1.5, col = 1, params = c(n_iter = 1e3))
  stereo_confidence(xpt$t, pch = 16, cex = 1.5, col = 2, params = c(n_iter = 1e3))
}))
```
