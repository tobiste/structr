# Centering vectors

Rotate vector object to position that ortensor eigenvectors are parallel
to axes of coordinate system: E3\|\|X (north-south), E2\|\|X(east-west),
E1\|\|X(vertical). Useful when one wants to inspect the distribution of
vectors, especially when vectors plot near the perimeter of the
stereonet.

## Usage

``` r
center(x, max_vertical = FALSE)
```

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

- max_vertical:

  Whether the maximum of the von Mises-Fisher distribution is already
  vertical or not.

## Value

Object of class of `x`

## See also

[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)

## Examples

``` r
set.seed(1)
mu <- Line(120, 10)
x <- rkent(100, mu = mu, k = 20, b = 5)
x_centered <- center(x)

# plot results
plot(x, col = "grey")
points(x_centered, col = "#B63679", pch = 16)
legend("topright", legend = c("original", "centered"), col = c("grey", "#B63679"), pch = 16)
```
