# Direct inversion For Fault slip in \\S^5\\

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

  6D unit y-vector on \\S^5\\ representing the tensor

- `alpha`:

  per-fault angular misfit (unsigned, 0-90°)

- `mean_alpha`:

  mean angular misfit across all faults

## Details

A reduced stress tensor is represented as a unit 5-vector \\x\\ on the
5-sphere \\S^5\\, embedded in the 6D space of symmetric tensors via the
\\y\\-vector (Eq. 6):

\$\$y = (1/\sqrt{2}) \times (s\_{11}, s\_{22}, s\_{33}, \sqrt{2} \times
s\_{23}, \sqrt{2} \times s\_{31}, \sqrt{2} \times s\_{12})\$\$

with normalisation constraints (Eqs. 11-12):

\\\sigma_I = s\_{11} + s\_{22} + s\_{33} = 0\\ (deviator)

\\\sigma\_{II} = (s\_{11}^2+s\_{22}^2+s\_{33}^2)/2 +
s\_{23}^2+s\_{31}^2+s\_{12}^2 = 1\\ (unit norm)

These constraints confine \\y\\ to the intersection of a hyperplane and
the 6-sphere, which defines the 5-sphere \\S^5\\.

For each fault datum \\\\n, \text{slip}\\\\, define the 6-vector
\\\epsilon'\\ (Eq. 8, derived from Eq. 4) such that the Wallace-Bott
condition becomes:

\\\epsilon'\_i \cdot y = 0\\ (shear stress perpendicular to \\b_i = n_i
\times \text{slip}\_i\\)

The optimal \\y\\ minimises \\\sum_i (\epsilon'\_i \cdot y)^2\\, i.e.
the smallest-eigenvalue eigenvector of \\M = \sum_i \epsilon'\_i \otimes
\epsilon'\_i\\.

The computation is carried out in the 5D deviator subspace via an
explicit orthonormal basis \\B\\ (\\6 \times 5\\), constructed as the
complement of the deviator null direction \\h =
(1,1,1,0,0,0)/\sqrt{3}\\, to avoid degeneracy in the \\6\times 6\\
projected matrix.

The sense condition (Eq. 5 / Eq. 10) selects between the two antipodal
solutions via a majority vote on \\\tau_i \cdot \text{slip}\_i \> 0\\
(positive work).

### Distance metrics

*Angular stress distance* (or *Lode angle*) \\\Theta\\ (Eq. 19):
great-circle distance on \\S^5\\, \\\Theta = \arccos{\|y_1 \dot
y_2\|})\\, range \\\[0, 90\]\\ deg for reduced tensors.

*Orife-Lisle distance* \\D\\ (Eq. 26): Euclidean chord on \\S^5\\, \\D =
2 \sin(\Theta/2)\\, range \\\[0, 2\]\\.

*Michael distance* \\D_M\\ (Eq. 34): \\D_M = 1 - \cos(\Theta)\\, range
\\\[0, 2\]\\.

Relationship to misfit angles (Eq. 37): \\\bar{d}\\ approx \\\Theta\\,
where \\\bar{d}\\ is the mean angular misfit between predicted and
observed slip directions.

## References

Yamaji, A., & Sato, K. (2006). Distances for the solutions of stress
tensor inversion in relation to misfit angles that accompany the
solutions. Geophysical Journal International, 167(2), 933–942.
[doi:10.1111/j.1365-246X.2006.03188.x](https://doi.org/10.1111/j.1365-246X.2006.03188.x)

## See also

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md),
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md)

Other yamaji:
[`slip_inversion_yamaji_sato_boot()`](https://tobiste.github.io/structr/reference/slip_inversion_yamaji_sato_boot.md)

## Examples

``` r
set.seed(20250411)

nx <- length(angelier1990)
par(mfrow = c(2, nx / 2))

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

dev.off()
#> null device 
#>           1 
```
