# Rigid Grain Net (RGN)

The rigid grain net (Jessup et al. 2007) plots the distribution of the
finite strain ratio \\R\\ and orientation \\\theta\\ of porphyroclasts
over the theoretical distribution of tailless clasts. The plot estimates
the critical shape factor \\B^\ast\\ marking the transition between
stable-sink positions and infinitely rotating porphyroclasts. This
critical shape factor is interpreted as the mean kinematic vorticity
number \\W_m\\.

## Usage

``` r
RGN_plot(
  r,
  theta,
  angle_error = 3,
  n_iter = 2000L,
  probs = 0.975,
  p_min = 0.95,
  grid = 0.05,
  main = "Rigid-Grain-Net",
  guides = NULL,
  ...
)
```

## Arguments

- r:

  numeric. Porphyroclast aspect ratio (long axis / short axis).

- theta:

  numeric. Angle between long axis and foliation (degrees).

- angle_error:

  numeric. Standard deviation of the angle measurement error in degrees.
  `3` by default.

- n_iter:

  integer. Number of bootstrap resamples. `2000` by default.

- probs:

  numeric. Probability in \\\[0, 1\]\\ used to estimate the critical
  shape factor, i.e. the quantile of shape factors of clasts classified
  as outside the critical hyperbola. `0.975` by default.

- p_min:

  numeric. Minimum probability that a clast lies outside the critical
  threshold for it to be used in the estimate. `0.95` by default.

- grid:

  numeric. Spacing of hyperbolae.

- main:

  character. The main title (on top).

- guides:

  logical. Whether guides of the projection grid should be added to the
  plot. Defaults to `getOption("structr.guides")`.

- ...:

  plotting arguments passed to
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html)

## Value

A plot. Invisibly, a list with `values` (matrix of `B` and `theta`),
`Wm` (point estimate of the critical shape factor) and `Wm_CI`
(bootstrap percentile interval).

## Details

### The Mean Kinematic Vorticity Number

The mean kinematic vorticity number \\W_m\\ is estimated as the upper
`probs` quantile of the shape factor \\B^\ast\\ of clasts that lie
outside the critical hyperbola (Jessup et al. 2007, Eq. 1, evaluated at
\\W_m = B^\ast\\).

The orientation \\\theta\\ of the porphyroclast long axis relative to
the foliation is given by \$\$\theta = \frac{1}{2} \sin^{-1}\Bigl\\
\frac{W_m}{B^{\ast}} \left( \sqrt{1-W_m^2} - \sqrt{{B^{\ast}}^2 - W_m^2}
\right) \Bigr\\\$\$ where \\B^{\ast}\\ is the shape factor \$\$B^{\ast}
= \frac{M_x^2 - M_n^2}{M_x^2 + M_n^2}\$\$ and \\M_n\\ and \\M_x\\ are
the short and long axes of the clast (Jessup et al. 2007).

### Confidence Interval of the Estimate

Measurement error in \\\theta\\ is handled probabilistically rather than
with a fixed margin. For each clast, the probability that its true angle
lies outside the critical threshold is \$\$p_i =
\Phi\left(\frac{\|\theta_i\| -
\theta\_{c,i}}{\sigma\_\theta}\right),\$\$ where \\\theta\_{c,i}\\ is
the threshold at that clast's \\B^\ast\\ and \\\sigma\_\theta\\ is
`angle_error`. Only clasts with \\p_i \ge\\ `p_min` are kept. Near the
vertex of the hyperbola the threshold is small, so a few degrees of
error can move a clast across it. Clasts with large \\B^\ast\\ and
uncertain classification are excluded rather than allowed to set the
upper quantile.

The confidence interval is a percentile interval from a non-parametric
bootstrap over clasts (`n_iter` resamples; after Stephan et al (2025).
Measurement error enters through the classification probabilities only,
so it is not added again through random perturbation of the data.
Because the estimate depends on `p_min`, the sensitivity of the result
to this choice should be reported by repeating the estimate over a range
of values

## References

Jessup, M. J., Law, R. D., & Frassi, C. (2007). The rigid grain net
(RGN): an alternative method for estimating mean kinematic vorticity
number (\\W_m\\). *Journal of Structural Geology*, 29(3), 411–421.
[doi:10.1016/j.jsg.2006.11.003](https://doi.org/10.1016/j.jsg.2006.11.003)

Stephan, T., Phillips, N., Tiitto, H., Perez, A., Nwakanma, M., Creaser,
R., & Hollings, P. (2025). Going with the flow - Changes of vorticity
control gold enrichment in Archean shear zones (Shebandowan Greenstone
Belt, Superior Province, Canada). *Journal of Structural Geology*, 201,
105542.
[doi:10.1016/j.jsg.2025.105542](https://doi.org/10.1016/j.jsg.2025.105542)

## Examples

``` r
data(shebandowan)
set.seed(20250411)
RGN_plot(shebandowan$r, shebandowan$phi, col = "darkred", guides = TRUE, n_iter = 1e3)
```
