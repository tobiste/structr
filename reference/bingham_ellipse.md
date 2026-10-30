# Confidence ellipse of an Eigenvector of the orientation tensor

Confidence ellipse of an Eigenvector of the orientation tensor

## Usage

``` r
bingham_ellipse(lambda, nu, cross_moments, n, level = 0.95, which = 1L)
```

## Source

modified after the python package `apsg` (v.2.0.6) by O. Lexa

## Arguments

- lambda:

  numeric. Eigenvalues of the orientation tensor, descending. See
  [`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)
  for details

- nu:

  Corresponding Eigenvectors. Object of class `Spherical`. See
  [`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)
  for details

- cross_moments:

  3x3 matrix (see
  [`axial_cross_moments()`](https://tobiste.github.io/structr/reference/axial_cross_moments.md))

- n:

  integer. Sample size

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

- which:

  integer. Index of the Eigenvector of the ellipse is centered on. One
  of `1` (the default), `2` or `3`.

## Value

list

## Details

Large-sample (central limit theorem) approximation: the angle of the
sample eigenvector \\\nu_i\\ toward \\\nu_j\\ is approximately normal
with variance \\m\_{ij} / (n (\lambda_i - \lambda_j)^2)\\, with
\\m\_{ij} = \text{mean}((x \cdot \nu_i)^2 (x \cdot \nu_j)^2)\\. The
semi-angles of the region are \$\$\arcsin(\sqrt(\chi^2_2(\text{level})
\times m\_{ij} / (n (\lambda_i - \lambda_j)^2)))\$\$, 90° for equal
eigenvalues.

The moments make it valid for any distribution of axes; for rotationally
symmetric data it is the cone of FLE eq. 6.11-6.13. Under the Bingham
model \\m\_{ij} / (\lambda_i - \lambda_j)^2\\ tends to \\1 / (2
\Delta\_{ij})\\, the semi-axes \\ \sqrt{\big(\chi^2_2(\text{level}) / (2
n \Delta\_{ij})\big)}\\ of Bingham (1974, eq. 7.24).

## See also

Other confidence:
[`bingham_confidence()`](https://tobiste.github.io/structr/reference/bingham_confidence.md),
[`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md),
[`fisher_confidence()`](https://tobiste.github.io/structr/reference/fisher_confidence.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md)

## Examples

``` r
x <- example_planes
t <- ot_eigen(x)
m <- axial_cross_moments(x, ot_eigen(x)$vectors)
r4 <- bingham_ellipse(lambda = t$values, nu = t$vectors, cross_moments = m, n = nrow(x), which = 1)
print(r4)
#> $mu
#> Vector (Vec3) object (n = 1):
#>          x          y          z 
#> -0.9257945  0.1785999  0.3331765 
#> 
#> $level
#> [1] 0.95
#> 
#> $axes
#> Vector (Vec3) object (n = 2):
#>              x          y         z
#> [1,] 0.2407711 -0.4008846 0.8839235
#> [2,] 0.2914340  0.8985508 0.3281351
#> 
#> $gamma
#> [1] 9.028139 3.794581
#> 

stereoplot()
points(x, col = "grey", pch = 16)
stereo_confidence(bingham_confidence(r4), col = "red")
```
