# Stress tensor inversion via the Yamaji and Sato (2006) eigenvector method.

Stress tensor inversion via the Yamaji and Sato (2006) eigenvector
method.

## Usage

``` r
slip_inversion_yamaji_sato(x, weights = NULL, flip = FALSE)
```

## Arguments

- x:

  object of class `"Pair"` or `"Fault"` with at least 4 rows.

- weights:

  numeric. Weightings for the faults. Must have the same length as `x`

- flip:

  logical. Flip if you want to have the negative stress tensor, i.e.
  sigma 1 and 3 will be flipped.

## Value

Same output as
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md)
plus

- `y`:

  6D unit y-vector on S^5 representing the tensor

- `alpha`:

  per-fault angular misfit (unsigned, 0-90°)

- `mean_alpha`:

  mean angular misfit across all faults

## References

Yamaji, A., & Sato, K. (2006). Distances for the solutions of stress
tensor inversion in relation to misfit angles that accompany the
solutions. Geophysical Journal International, 167(2), 933–942.
https://doi.org/10.1111/j.1365-246X.2006.03188.x

## See also

[`slip_inversion_yamaji_sato_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_yamaji_sato_boot.md)

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_hansen_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen_boot.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md),
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md)

## Examples

``` r
set.seed(20250411)
nx <- length(angelier1990)
par(mfrow = c(1, nx))

invisible(lapply(seq_len(nx), function(i) {
  # inversion
  x <- angelier1990[[i]]
  res <- slip_inversion_yamaji_sato(x)

  # some stress shape
  phi_val <- round(res$stress_shape$phi, 2)

  # misfit
  rup_val <- round(res$misfit$rup_mean, 2)

  # Plot the faults (color-coded by RUP%) and show the principal stress axes
  stereoplot(title = names(angelier1990)[i], guides = FALSE)
  stereo_shmax(res$SHmax)
  fault_plot(x, col = assign_col(res$misfit$rup))
  points(res$principal_axes, col = 1:3, pch = 16, cex = 1.5)
  text(res$principal_axes,
    label = rownames(res$principal_axes),
    col = 1:3, adj = -.25
  )
  legend("topleft", col = 2:4, legend = rownames(res$principal_axes), pch = 16)
  title(sub = bquote(Phi == .(phi_val) ~ "|" ~ bar("RUP") == .(rup_val) * "%"))
}))
```
