# Confidence region for Fisher statistics

Circle of angular radius \\\alpha\\ (degrees) around the mean vector
\\\mu\\.

## Usage

``` r
fisher_confidence(mu, alpha, n = 512L)
```

## Arguments

- mu:

  Object of class `spherical`

- alpha:

  numeric. Angle in degrees.

- n:

  integer. Resolution of the cone (number of points for the cone line)

## Value

object of class `"confidence_region"`

## See also

Other confidence:
[`bingham_confidence()`](https://tobiste.github.io/structr/reference/bingham_confidence.md),
[`bingham_ellipse()`](https://tobiste.github.io/structr/reference/bingham_ellipse.md),
[`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md)

## Examples

``` r
set.seed(20250411)
x <- rvmf(n = 100, Line(120, 55), k = 10)
mu <- sph_mean(x)
alpha <- sph_confidence_angle(x)

r1 <- fisher_confidence(mu, alpha = alpha)

plot(x, col = "grey")
points(mu, col = "red", pch = 16, cex = 1.5)
stereo_confidence(r1, col = "red")
```
