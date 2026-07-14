# Stress Inversion for Fault-Slip Data

Convenience function for direction inversion of fault-slip data to
derive the reduced stress tensor.

## Usage

``` r
slip_inversion(
  x,
  method = c("michael", "angelier", "hansen", "yamaji", "wissi"),
  ...
)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates. Object must be complete, i.e. no `NA` values. For
  Michael's, Angelier's, and Yamaji-Sato's methods, at least 4 rows of
  fault measurements are required, while Hansen's method requires at
  least 7.

- method:

  character. The inversion algorithm, one of `"michael"` (the default)
  for a bootstrapped linear inversion after Micheal (1984), `"angelier"`
  for an iterative direct inversion after Angelier (1990) and Mostafa
  (2005), `"yamaji"` for direct inversion using the 5d parameter space
  after Yamaji and Sato (2006), and `"hansen"` for direct inversion
  using the 9d parameter space after Hansen (2013).

- ...:

  arguments passed to
  [`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
  [`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
  or
  [`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md)
  depending on `method`.

## Value

a named list with the following components:

- `stress_tensor`:

  `"ellipsoid"` object. Best-fit devitoric stress tensor in input
  coordinate frame

- `principal_axes`:

  `"Line"` objects. Orientation of the principal stress axes as unit
  vectors (max to min)

- `tensor_params`:

  the four tensor parameters (Eq. 4.87)

- `principal_vals`:

  eigenvalues of the stress tensor (\\\sigma_1 \>= \sigma_2 \>=
  \sigma_3\\)

- `stress_shape`:

  list Stress shape ratio. See
  [`stress_shape()`](https://tobiste.github.io/structr/reference/stress_shape.md).

- `misfit`:

  list. Misfit parameters. See
  [`slip_inversion_misfit()`](https://tobiste.github.io/structr/reference/slip_inversion_misfit.md).

- `SHmax`:

  numeric. Direction of maximum horizontal stress (in degrees)

- `tau_mean`:

  numeric. Average resolved shear stress on each plane. Should be close
  to 1.

- `stress_components`:

  matrix. The resolved shear and normal stresses, the slip and dilation
  tendency on each plane. See
  [`tau2shearnorm()`](https://tobiste.github.io/structr/reference/tau-comp.md)
  and
  [`tau2tendency()`](https://tobiste.github.io/structr/reference/tau-comp.md).

- `n_iter`:

  number of Mostafa iterations performed

- `method`:

  character. The inversion method used, equal to `method` argument.

## References

Angelier, J. (1990). Inversion of field data in fault tectonics to
obtain the regional stress—III. A new rapid direct inversion method by
analytical means. Geophys. J. Int, 103, 363–376.
<https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

Hansen, J. A. (2013). Direct inversion of stress, strain or strain rate
including vorticity: A linear method of homogenous fault-slip data
inversion independent of adopted hypothesis. Journal of Structural
Geology, 51, 3–13. https://doi.org/10.1016/j.jsg.2013.03.014

Michael, A. J. (1984). Determination of stress from slip data: Faults
and folds. Journal of Geophysical Research: Solid Earth, 89(B13),
11517–11526. <https://doi.org/10.1029/JB089iB13p11517>

Yamaji, A., & Sato, K. (2006). Distances for the solutions of stress
tensor inversion in relation to misfit angles that accompany the
solutions. Geophysical Journal International, 167(2), 933–942.
https://doi.org/10.1111/j.1365-246X.2006.03188.x

## See also

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_hansen_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen_boot.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md),
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md),
[`slip_inversion_yamaji_sato()`](https://tobiste.github.io/structr/reference/slip_inversion_yamaji_sato.md)

## Examples

``` r
set.seed(20250411)
# Use Angelier examples
par(mfrow = c(1, length(angelier1990)))
invisible(lapply(angelier1990, function(x) {
  # Inversion after Michael (1984)
  res_michael <- slip_inversion(x, method = "michael", n_iter = 100, n = 100, res = 100)

  # Inversion after Angelier (1990)
  res_angelier <- slip_inversion(x, method = "angelier")
  
  res_yamaji <- slip_inversion(x, method = "yamaji")

  res_hansen <- slip_inversion(x, method = "hansen", type = "6d")
  
  res_wissi <- slip_inversion(x, method = 'wissi')

  stereoplot(guides = FALSE)
  fault_plot(x, col = "gray80")
  points(res_michael$principal_axes, pch = 1:3, col = 2)
  points(res_angelier$principal_axes, pch = 1:3, col = 3)
  points(res_yamaji$principal_axes, pch = 1:3, col = 4)
  points(res_hansen$principal_axes, pch = 1:3, col = 5)
  points(res_wissi$principal_axes, pch = 1:3, col = 6)
  legend("topleft",
    pch = 1,
    legend = c("Michael (1984)", "Angelier (1990)", "Yamaji & Sato (2006)", "Hansen (2013)", "WISSI"),
    col = 2:6
  )
  legend("bottomright",
    pch = 1:3,
    legend = c("S1", "S2", "S3")
  )
}))
#> Error in normals %*% TR: requires numeric/complex matrix/vector arguments
```
