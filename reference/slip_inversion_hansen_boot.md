# 9D Direct Inversion for Fault Slip Including Vorticity with Confidence Intervals

Fault-slip inversion method after Hansen (2013) with bootstrapped
confidence intervals

## Usage

``` r
slip_inversion_hansen_boot(
  x,
  friction = 0.6,
  flip = FALSE,
  type = c("9d", "6d"),
  n_iter = 100L,
  conf.level = 0.95,
  ...
)
```

## Arguments

- x:

  object of class `"Pair"` or `"Fault"` with at least 7 rows.

- friction:

  numeric. Coefficient of friction (0.6 by default)

- flip:

  logical. Flip if you want to have the negative stress tensor, i.e.
  sigma 1 and 3 will be flipped.

- type:

  character. Inversion method, either `"9d"` (the default) for using the
  9-dimensional or `"6d"` for the 6-dimensional parameter space.

- n_iter:

  integer. Number of bootstrap replicates (100 by default)

- conf.level:

  numeric. Confidence level of the interval (0.95 by default)

- ...:

  optional parameters passed to
  [`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md)

## Value

See
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md)
and
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md)

## See also

Other stress-inversion:
[`Fault_PT()`](https://tobiste.github.io/structr/reference/Fault_PT.md),
[`slip_inversion()`](https://tobiste.github.io/structr/reference/slip_inversion.md),
[`slip_inversion_angelier()`](https://tobiste.github.io/structr/reference/slip_inversion_angelier.md),
[`slip_inversion_hansen()`](https://tobiste.github.io/structr/reference/slip_inversion_hansen.md),
[`slip_inversion_michael()`](https://tobiste.github.io/structr/reference/slip_inversion_michael.md),
[`slip_inversion_simple()`](https://tobiste.github.io/structr/reference/slip_inversion_simple.md),
[`slip_inversion_wissi()`](https://tobiste.github.io/structr/reference/slip_inversion_wissi.md),
[`slip_inversion_yamaji_sato()`](https://tobiste.github.io/structr/reference/slip_inversion_yamaji_sato.md)

## Examples

``` r
set.seed(20250411)
res <- slip_inversion_hansen_boot(osmundsen2010, n_iter = 100, n = 1000, res = 100)

# some stress shape
phi_val <- round(res$phi_CI, 2)

# vorticity
w_val <- round(res$vorticity_mag_CI, 2)

# Plot the faults
plot(osmundsen2010, col = "grey", lwd = 0.1, cex = 0.5)
stereo_confidence(res$principal_axes_CI$sigma1, col = 2)
stereo_confidence(res$principal_axes_CI$sigma2, col = 3)
stereo_confidence(res$principal_axes_CI$sigma3, col = 4)
stereo_confidence(res$vorticity_axis_CI, col = 5)
text(res$principal_axes, label = rownames(res$principal_axes), col = 2:4, adj = -.25)
text(res$vorticity_axis, labels = bquote(omega), col = 5, adj = -.5)
title(
  main = "Lofoten / Northern Norway\n(Osmundsen et al. 2010)",
  sub = bquote(atop(
    varphi ~ "(95% CI)" == "[" * .(phi_val[1]) * "," ~ .(phi_val[2]) * "]",
    ~omega ~ "(95%)" == "[" * .(w_val[1]) * "," ~ .(w_val[2]) * "]"
  ))
)
```
