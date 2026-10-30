# Confidence region for Bingham statistics

Bingham confidence ellipse

## Usage

``` r
bingham_confidence(stats, n = 512)
```

## Arguments

- stats:

  Output of
  [`bingham_statistics()`](https://tobiste.github.io/structr/reference/bingham_statistics.md)

- n:

  integer. Resolution of the cone (number of points for the cone line)

## Value

object of class `"confidence_region"`

## See also

Other confidence:
[`bingham_ellipse()`](https://tobiste.github.io/structr/reference/bingham_ellipse.md),
[`confidence_ellipse()`](https://tobiste.github.io/structr/reference/confidence_ellipse.md),
[`fisher_confidence()`](https://tobiste.github.io/structr/reference/fisher_confidence.md),
[`stereo_confidence()`](https://tobiste.github.io/structr/reference/stereo_confidence.md)

## Examples

``` r
set.seed(20250411)
a <- cov(iris[, 1:3])
x <- rbingham(100, a, "Line")
r <- bingham_statistics(x)
r2 <- bingham_confidence(r)
plot(x, col = "grey")
points(r$mu, col = "red", pch = 16, cex = 1.5)
stereo_lines(r2$ellipse, col = "red")


r3 <- bingham_statistics(example_planes)
stereoplot()
points(example_planes, col = "grey", pch = 16)
stereo_confidence(bingham_confidence(r3), col = "red")
```
