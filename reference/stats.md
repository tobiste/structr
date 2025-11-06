# Statistical estimators of the distribution of a set of vectors

Statistical estimators of the distribution of a set of vectors

## Usage

``` r
sph_mean(x, na.rm = TRUE, ...)

sph_sd(x, ...)

sph_var(x, ...)

sph_confidence_angle(x, w = NULL, alpha = 0.05, na.rm = TRUE)

rdegree(x, w = NULL, na.rm = FALSE)

sd_error(x, w = NULL, na.rm = FALSE)

delta(x, w = NULL, na.rm = TRUE)

estimate_k(x, w = NULL, na.rm = FALSE)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- na.rm:

  logical. Whether `NA` values should be removed before the computation
  proceeds.

- ...:

  arguments passed to function call

- w:

  numeric. Optional weights for each observation.

- alpha:

  numeric. Significance level for the confidence angle (default is 0.05
  for a 95% confidence angle).

## Details

These statistical estimators are based on the resultant vector of a set
of \\n\\ vectors \\x_1, \ldots, x_n\\ (Mardia 1972). The resultant
vector is given by \$\$\bar{\mathbf{x}} = \sum\_{i=1}^{n}
\mathbf{x}\_i\$\$

The mean resultant is defined as \$\$\bar{\mathbf{R}} =
\|\|\bar{\mathbf{x}}\|\| = \sqrt{x_x^2 + x_y^2 + x_z^2}\$\$

`sph_mean` returns the spherical mean of a set of vectors (object of
class of `x`), that is the maximum likelihood estimate of the mean
direction given by the arithmetic mean of the vector components
normalized by the mean resultant vector: \$\$\mu =
\frac{\bar{\mathbf{x}}}{\bar{\mathbf{R}}}\$\$

`sph_var` returns the spherical variance (numeric). \$\$S = 1 -
\bar{\mathbf{R}}\$\$

`sph_sd` returns the spherical standard deviation (numeric) given as the
half apical angle of a cone about the mean vector. In degrees if `x` is
a `"Plane"` or `"Line"`, or in radians if otherwise. \$\$s =
\sqrt{\log(1 / \bar{\mathbf{R}}^2))}\$\$

`delta` returns the half apical angle of the cone containing ~63% of the
data (in degrees if `x` is a `"Plane"` or `"Line"`, or in radians if
otherwise). For enough large sample it approaches the angular standard
deviation (`"csd"`) of the Fisher statistics. \$\$\delta =
\arccos(\bar{\mathbf{R}})\$\$

`rdegree` returns the degree of preferred orientation of vectors, range:
(0, 1). \$\$r = \frac{2 \bar{\mathbf{R}} - n}{n}\$\$

`sd_error` returns the spherical standard error (numeric). If the number
of data is less than 25, if will print a additional message, that the
output value might not be a good estimator. \$\$\text{SDE} =
\sqrt{\frac{1 - \frac{1}{n} \sum\_{i=1}^{n} (\mu \cdot x_i)^2}{n
\bar{\mathbf{R}}^2}}\$\$

`sph_confidence_angle` returns the half-apical angle \\q\\ of a cone
about the mean \\\mu\\ (in degrees if `x` is a `"Plane"` or `"Line"`, or
in radians if otherwise). The \\100(1-\alpha)\\\\ confidence interval is
than given by \\\mu \pm q\\. \$\$q = \arcsin(\sqrt{-\log(\alpha)} \cdot
\text{SDE})\$\$

`estimate_k` returns the estimated concentration of the von Mises-Fisher
distribution \\\hat{\kappa}\\ (after Sra, 2011). \$\$\hat{\kappa} =
\frac{\bar{R}(p - \bar{R}^2)}{1 - \bar{R}^2}\$\$ where \\p\\ is the
dimension of the data (3 for spherical data).

## References

Mardia, Kanti; Jupp, P. E. (1999). Directional Statistics. John Wiley &
Sons Ltd. ISBN 978-0-471-95333-3.

Sra, S. A short note on parameter approximation for von Mises-Fisher
distributions: and a fast implementation of I s (x). Comput Stat 27,
177–190 (2012). https://doi.org/10.1007/s00180-011-0232-x

## See also

[`projected_mean()`](https://tobiste.github.io/structr/reference/projected_mean.md)
for projected mean,
[`geodesic_mean_line()`](https://tobiste.github.io/structr/reference/geodesic-line.md)
and
[`geodesic_mean_pair()`](https://tobiste.github.io/structr/reference/mean-pair.md)
for the geodesic mean of lines and pairs, respectively. or
[`geodesic_mean()`](https://tobiste.github.io/structr/reference/geodesic-mean.md)
and
[`geodesic_var()`](https://tobiste.github.io/structr/reference/geodesic-var.md)
as a convenience wrapper for all spherical data types.

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Line(120, 50), k = 5)
sph_mean(x)
#> Line object (n = 1):
#>   azimuth    plunge 
#> 115.21007  52.67744 
sph_sd(x)
#> [1] 40.82336
sph_var(x)
#> [1] 0.224176
delta(x)
#> [1] 39.1202
rdegree(x)
#> [1] -0.9844835
sd_error(x)
#> [1] 4.01451
sph_confidence_angle(x)
#> [1] 6.965534
estimate_k(x)
#> [1] 4.673486

#' weights:
x2 <- Line(c(0, 0), c(0, 90))
sph_mean(x2)
#> Line object (n = 1):
#> azimuth  plunge 
#>       0      45 
sph_mean(x2, w = c(1, 2))
#> Line object (n = 1):
#>  azimuth   plunge 
#>  0.00000 63.43495 
sph_var(x2)
#> [1] 0.2928932
sph_var(x2, w = c(1, 2))
#> [1] 0.254644
```
