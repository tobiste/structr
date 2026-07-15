# Weighted Iterative Sigma-Space Inversion (WISSI)

Combines the four classic fault slip inversion algorithms (Michael,
1983; Angelier, 1990; Yamaji and Sato, 2006; and Hansen, 2013) into a
single coherent framework operating in the Yamaji-Sato 5-sphere
sigma-space.

## Usage

``` r
slip_inversion_wissi(
  x,
  weights = NULL,
  sigma_alpha_deg = 10,
  gamma_max = 10,
  n_anneal = 8L,
  max_iter = 50L,
  tol_deg = 1e-04,
  run_stage4 = TRUE
)
```

## Arguments

- x:

  object of class `"Pair"` or `"Fault"` with at least 4 rows.

- weights:

  numeric. Weightings for the faults. Must have the same length as `x`

- sigma_alpha_deg:

  Estimated slip direction measurement error in degrees. Used in Stage 4
  analytic uncertainty. Default `10.`

- gamma_max:

  Maximum sense annealing sharpness parameter. Higher values commit more
  strongly to the predicted slip sense. Default `10.` Set to `0` to
  disable sense annealing (fully sense-agnostic, like Hansen 2013).

- n_anneal:

  Number of annealing steps (outer loop of Stage 3). Default `8`.

- max_iter:

  Maximum inner iterations per annealing step. Default `50.`

- tol_deg:

  Convergence tolerance in angular stress distance (degrees). Default
  `1e-4`.

- run_stage4:

  Logical. Compute analytic uncertainty (Stage 4). Default `TRUE.`

## Value

A named list with:

- `stress_tensor`:

  3x3 reduced stress tensor (Cartesian frame)

- `y`:

  6D unit y-vector on \\S^5\\

- `principal_axes`:

  unit vectors of principal stress axes (max to min)

- `principal_vals`:

  eigenvalues of `stress_tensor` (decreasing)

- `alpha`:

  per-fault angular misfit (unsigned, 0-90°)

- `alpha_signed`:

  per-fault signed misfit (0-180°)

- `mean_alpha`:

  mean angular misfit across all faults (°)

- `suspected_flipped`:

  row indices where `alpha_signed` \> 90°

- `n_flipped_sense`:

  number of faults whose sense was corrected in Stage 3

- `slips_corrected`:

  sense-corrected slip matrix used in final inversion

- `mu`:

  per-fault magnitude weights from Stage 2/3

- `phi_sense`:

  per-fault tanh sense confidence from Stage 3

- `eigenvalue_gap`:

  \\\lambda_2 - \lambda_1\\ of \\M5\\ (condition number proxy)

- `M5_eigvals`:

  all 5 eigenvalues of final \\M5\\

- `unc`:

  Stage 4 uncertainty list (if `run_stage4 = TRUE`): `Cov5`: 5x5
  covariance matrix in sigma-space; `Cov_y6`: 6x6 covariance matrix
  (`y`-space); `eigval_gap`: eigenvalue gap (same as above);
  `cov_eigvals`: eigenvalues of `Cov5`; `sigma1_unc`: approx 1\\\sigma\\
  uncertainty on \\\sigma_1\\ orientation `Phi_unc`: approx 1\\sigma\\
  uncertainty on \\\phi\\

- `n_iter_total`:

  total number of inner iterations

## References

Stephan (in prep.)

## See also

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_hansen_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen_boot.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md),
[`slip_inversion_yamaji_sato()`](https://tobiste.github.io/structr/reference/slip_inversion_yamaji_sato.md)

## Examples

``` r
set.seed(20250411)

nx <- length(angelier1990)
par(mfrow = c(2, nx/2))

invisible(lapply(seq_len(nx), function(i) {
  # inversion
  x <- angelier1990[[i]]
  res <- slip_inversion_wissi(x)

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
