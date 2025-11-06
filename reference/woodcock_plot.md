# Fabric plot of Woodcock (1977)

Creates a fabric plot using the eigenvalue method

## Usage

``` r
woodcock_plot(
  x,
  labels = NULL,
  add = FALSE,
  max = 7,
  main = "Woodcock diagram",
  ...
)
```

## Arguments

- x:

  either an object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`,
  `"Pair"`, or `"Fault"` where the rows are the observations and the
  columns are the coordinates, or an `"ortensor"` object.

- labels:

  character. text labels

- add:

  logical. Should data be plotted to an existing plot?

- max:

  numeric. Maximum value for x and y axes. If `NULL`, it is calculated
  from the data.

- main:

  character. The main title for the plot.

- ...:

  optional graphical parameters

## Value

A plot and when stored as an object, the orientation tensor's
eigenvalues and eigenvectors as a list.

## References

Woodcock, N. H. (1977). Specification of fabric shapes using an
eigenvalue method. Geological Society of America Bulletin88,
1231\<U+2013\>1236.
http://pubs.geoscienceworld.org/gsa/gsabulletin/article-pdf/88/9/1231/3418366/i0016-7606-88-9-1231.pdf

## See also

[`vollmer()`](https://tobiste.github.io/structr/reference/vollmer.md),
[`ot_eigen()`](https://tobiste.github.io/structr/reference/ot_eigen.md)

Other fabric-plot:
[`flinn_plot()`](https://tobiste.github.io/structr/reference/flinn_plot.md),
[`hsu_plot()`](https://tobiste.github.io/structr/reference/hsu_plot.md),
[`vollmer-plot`](https://tobiste.github.io/structr/reference/vollmer-plot.md)

## Examples

``` r
set.seed(20250411)
mu <- Line(120, 50)
x <- rvmf(100, mu = mu, k = 1)
woodcock_plot(x, lab = "x")
y <- rvmf(100, mu = mu, k = 20)
woodcock_plot(y, lab = "y", add = TRUE, col = "red")
```
