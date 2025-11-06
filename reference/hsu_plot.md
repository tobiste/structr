# Fabric plot of Hsu (1965)

3D strain diagram using the Hsu (1965) method to display the natural
octahedral strain (Nádai, 1950) and Lode's parameter (Lode, 1926).

## Usage

``` r
hsu_plot(x, ...)

# S3 method for class 'ortensor'
hsu_plot(x, labels = NULL, add = FALSE, es.max = 3, main = "Hsu diagram", ...)

# S3 method for class 'spherical'
hsu_plot(x, ...)

# S3 method for class 'ellispoid'
hsu_plot(x, ...)

# Default S3 method
hsu_plot(x, labels = NULL, add = FALSE, es.max = 3, main = "Hsu diagram", ...)

# S3 method for class 'list'
hsu_plot(x, labels = NULL, add = FALSE, es.max = 3, main = "Hsu diagram", ...)
```

## Arguments

- x:

  accepts the following objects: a two-column matrix where first column
  is the ratio of maximum strain and intermediate strain (X/Y) and
  second column is the the ratio of intermediate strain and minimum
  strain (Y/Z); objects of class `"Vec3"`, `"Line"`, `"Ray"`, or
  `"Plane"`; or `"ortensor"` objects.

- ...:

  plotting arguments passed to
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html)

- labels:

  character. text labels

- add:

  logical. Should data be plotted to an existing plot?

- es.max:

  maximum strain for scaling.

- main:

  character. The main title (on top).

## Value

plot and when stored as object, a list containing the Lode parameter
`lode` and the natural octahedral strain `es`.

## References

Lode, W. (1926). Versuche über den Einfluß der mittleren Hauptspannung
auf das Fließen der Metalle Eisen. Kupfer und Nickel. Zeitschrift Für
Physik, 36(11–12), 913–939.
[doi:10.1007/BF01400222](https://doi.org/10.1007/BF01400222)

Nádai, A. (1950). Theory of flow and fracture of solids. McGraw-Hill.

Hsu, T. C. (1966). The characteristics of coaxial and non-coaxial strain
paths. Journal of Strain Analysis, 1(3), 216–222.
[doi:10.1243/03093247V013216](https://doi.org/10.1243/03093247V013216)

Hossack, J. R. (1968). Pebble deformation and thrusting in the Bygdin
area (Southern Norway). Tectonophysics, 5(4), 315–339.
[doi:10.1016/0040-1951(68)90035-8](https://doi.org/10.1016/0040-1951%2868%2990035-8)

## See also

[`lode()`](https://tobiste.github.io/structr/reference/ellipsoid-params.md)
for Lode parameter, and
[nadai](https://tobiste.github.io/structr/reference/ellipsoid-params.md)
for natural octahedral strain.

Other fabric-plot:
[`flinn_plot()`](https://tobiste.github.io/structr/reference/flinn_plot.md),
[`vollmer-plot`](https://tobiste.github.io/structr/reference/vollmer-plot.md),
[`woodcock_plot()`](https://tobiste.github.io/structr/reference/woodcock_plot.md)

## Examples

``` r
R_XY <- holst[, "R_XY"]
R_YZ <- holst[, "R_YZ"]
hsu_plot(cbind(R_XY, R_YZ), col = "#B63679", pch = 16, type = "b")


set.seed(20250411)
mu <- Line(120, 50)
x <- rvmf(100, mu = mu, k = 1)
hsu_plot(x, labels = "x")

set.seed(20250411)
y <- rvmf(100, mu = mu, k = 20)
hsu_plot(ortensor(y), labels = "y", col = "red", add = TRUE)
```
