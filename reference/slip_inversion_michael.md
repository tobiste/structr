# Stress Inversion for Fault-Slip Data after Michael (1984)

Direct stress inversion (based on Michael, 1984) determines the
orientation of the principal stresses from fault slip data. Confidence
intervals are estimated by bootstrapping. This inversion is simplified
by the assumption that the magnitude of the tangential traction on the
various fault planes, at the time of rupture, is similar.

## Usage

``` r
slip_inversion_michael(
  x,
  n_iter = 100L,
  conf.level = 0.95,
  friction = 0.6,
  flip = FALSE,
  ...
)
```

## Arguments

- x:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

- n_iter:

  integer. Number of bootstrap samples (10 by default)

- conf.level:

  numeric. Confidence level of the interval (0.95 by default)

- friction:

  numeric. Coefficient of friction (0.6 by default)

- flip:

  logical. Flip if you want to have the negative stress tensor, i.e.
  sigma 1 and 3 will be flipped.

- ...:

  optional parameters passed to
  [`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)

## Value

Additionally, this child functions appends the following list
components:

- `principal_axes_CI`:

  list containing the confidesnce ellipses for the 3 principal stress
  vectors. See
  [`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)
  for details.

- `principal_vals_CI`:

  3-column vector containing the lower and upper margins of the
  confidence interval of the principal vals

- `SHmax_CI`:

  numeric. Confidence interval of `SHmax` angle

- `R_CI`,`phi_CI`,`bott_CI`:

  Confidence interval for `R`

- `alpha_CI`,`beta_CI`,`theta_CI`:

  numeric. Confidence intervals of `alpha`, `beta`, and `theta` angles

## Details

The goal of slip inversion is to find the single uniform stress tensor
that most likely caused the faulting events. With only slip data to
constrain the stress tensor the isotropic component can not be
determined, unless assumptions about the fracture criterion are made.
Hence inversion will be for the deviatoric stress tensor only. A single
fault can not completely constrain the deviatoric stress tensor a,
therefore it is necessary to simultaneously solve for a number of
faults, so that a single a that best satisfies all of the faults is
found.

## References

Michael, A. J. (1984). Determination of stress from slip data: Faults
and folds. Journal of Geophysical Research: Solid Earth, 89(B13),
11517–11526.
[doi:10.1029/JB089iB13p11517](https://doi.org/10.1029/JB089iB13p11517)

## See also

[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md)
for a simple P-T stress analysis,
[`SH()`](https://tobiste.github.io/structr/reference/SH.md) and
[`SH_from_tensor()`](https://tobiste.github.io/structr/reference/SH_from_tensor.md)
to calculate the azimuth of the maximum horizontal stress;
[`Mohr_plot()`](https://tobiste.github.io/structr/reference/Mohr_plot.md)
for graphical representation of the deviatoric stress tensor.

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md)

## Examples

``` r
# Use Angelier examples:
nx <- length(angelier1990)
par(mfrow = c(1, length(angelier1990)))

invisible(lapply(seq_len(nx), function(i){

# inversion
x <- angelier1990[[i]]
res <- slip_inversion_michael(x, n_iter = 100, n = 1000, res = 100)

# some stress shape
phi_val <- round(res$stress_shape$phi, 2)

# misfit
rup_val <- round(res$misfit$rup_mean, 2)

# Plot the faults (color-coded by RUP%) and show the principal stress axes
stereoplot(title = names(angelier1990)[i], guides = FALSE)
stereo_shmax(res$SHmax)
fault_plot(x, col = assign_col(res$misfit$rup))
stereo_confidence(res$principal_axes_CI$sigma1, col = 2)
stereo_confidence(res$principal_axes_CI$sigma2, col = 3)
stereo_confidence(res$principal_axes_CI$sigma3, col = 4)
text(res$principal_axes, label = rownames(res$principal_axes), col = 2:4, adj = -.25)
legend("topleft", col = 2:4, legend = rownames(res$principal_axes), pch = 16)
title(sub = bquote(varphi == .(phi_val) ~ "|" ~ bar("RUP") == .(rup_val) * '%'))
}))
```
