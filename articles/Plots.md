# Orientation Plots

This tutorial shows how to create various orientation plots using the
{structr} package, including stereographic and equal-area projections,
fabric plots, density plots, and fault plots.

``` r
library(structr)
```

Import and convert to spherical objects:

``` r
data(example_planes_df)
data(example_lines_df)

planes <- Plane(example_planes_df$dipdir, example_planes_df$dip)
lines <- Line(example_lines_df$trend, example_lines_df$plunge)
```

## Equal-area projection

Lambert equal area, lower hemisphere projection is the default plotting
setting.

``` r
stereoplot()
points(lines, col = "#B63679", pch = 19, cex = .5)
points(planes, col = "#000004", pch = 1, cex = .5)

legend("topright", legend = c("Lines", "Planes"), col = c("#B63679", "#000004"), pch = c(19, 1), cex = 1)
title(main = "Example data", sub = "Lambert equal area, lower hemisphere projection")
```

![Diagram showing the equal-area Lambert projection of some example
data](Plots_files/figure-html/stereo_earea-1.png)

> Points can be added using the
> [`points()`](https://rdrr.io/r/graphics/points.html) function.

## Stereographic projection

To change to equal angle stereographic, upper hemisphere projection,
just set the `earea` argument to `FALSE`, and the `upper.hem` argument
to `TRUE`:

``` r
stereoplot(earea = FALSE)
points(lines, col = "#B63679", pch = 19, cex = .5, earea = FALSE, upper.hem = TRUE)
points(planes, col = "#000004", pch = 1, cex = .5, earea = FALSE, upper.hem = TRUE)

legend("topright", legend = c("Lines", "Planes"), col = c("#B63679", "#000004"), pch = c(19, 1), cex = 1)
title(main = "Example data", sub = "Equal angle stereographic, upper hemisphere projection")
```

![Diagram showing the equal-angle stereographic projection of some
example data](Plots_files/figure-html/stereo_eangle-1.png)

### Great and small circles

Great and small circles can be added using the
[`lines()`](https://rdrr.io/r/graphics/lines.html) function.

Adding great circles for the first 10 vectors in planes:

``` r
stereoplot(guides = FALSE) # turn of guides for better visibility
lines(planes[1:10, ], col = "lightgrey", lty = 1)
points(planes[1:10, ], col = "#000004", pch = 1, cex = .5)
```

![Diagram showing
greatcircles](Plots_files/figure-html/stereo_greatcircle-1.png)

To plot a small circle with, e.g., a 10° radius, you need to specify the
`ang` argument in [`lines()`](https://rdrr.io/r/graphics/lines.html):

``` r
stereoplot(guides = FALSE)
points(lines[1:5, ], col = "#B63679", pch = 19, cex = .5)
lines(lines[1:5, ], ang = 10, col = "#B63679")
```

![Diagram showing
smallcircles](Plots_files/figure-html/stereo_smallcircle-1.png)

## Fabric plots

The **Eigenvalues** of the orientation tensor describe the shape of the
distribution of these vectors, that is who clustered, cylindrical or
random these vectors are distributed.

A Fabric plot visualizes the shape of the distribution by plotting the
eigenvalues of the orientation tensor. Three different diagram are
provided by {structr}, namely the triangular *Vollmer plot*[¹](#fn1),
the logarithmic biplot (*Woodcock plot*)[²](#fn2), and the Lode
parameter vs. natural octahedral strain diagram (*Hsu plot*)\[\|^hsu\].

### Vollmer plot

[`vollmer_plot()`](https://tobiste.github.io/structr/reference/vollmer-plot.md)
creates a triangular plot showing the shape of the orientation
distribution (after Vollmer, 1990).

``` r
vollmer_plot(planes, col = "#000004", pch = 1, cex = 2)
vollmer_plot(lines, add = TRUE, col = "#B63679", pch = 19, cex = 2)
legend("topright", legend = c("Lines", "Planes"), col = c("#B63679", "#000004"), pch = c(19, 1), cex = 1)
```

![Diagram showing eigen-vector analysis of some example data plotted in
an triangular plot after Vollmer
1990](Plots_files/figure-html/vollmer-1.png)

### Woodcock plot

[`woodcock_plot()`](https://tobiste.github.io/structr/reference/woodcock_plot.md)
creates a logarithmic biplot showing the shape of the orientation
distribution (after Woodcock, 1977).

``` r
woodcock_plot(planes, col = "#000004", pch = 1, cex = 2)
woodcock_plot(lines, add = TRUE, col = "#B63679", pch = 19, cex = 2)
legend("topright", legend = c("Lines", "Planes"), col = c("#B63679", "#000004"), pch = c(19, 1), cex = 1)
```

![Diagram showing eigen-vector analysis of some example data plotted in
an biplot plot after Woodcok
1977](Plots_files/figure-html/woodcock-1.png)

### Hsu plot

[`hsu_plot()`](https://tobiste.github.io/structr/reference/hsu_plot.md)
creates a Lode parameter[³](#fn3) vs. natural octahedral strain[⁴](#fn4)
diagram showing the shape of the orientation distribution (after Hsu,
1965).

``` r
hsu_plot(planes, col = "#000004", pch = 1, cex = 2)
hsu_plot(lines, add = TRUE, col = "#B63679", pch = 19, cex = 2)
legend("topright", legend = c("Lines", "Planes"), col = c("#B63679", "#000004"), pch = c(19, 1), cex = 1)
```

![Diagram showing eigen-vector analysis of some example data plotted in
a Hsu plot](Plots_files/figure-html/hsu-1.png)

## Density plots

**Kamb** contours[⁵](#fn5) and densities can be added to an existing
projection plot using the `contour` functions. Weighted densities can be
controlled by the `weights` argument and are useful when the orientation
measurements have different accuracies.

``` r
example_planes_df$quality <- ifelse(is.na(example_planes_df$quality), 6, example_planes_df$quality) # replacing NA values with 6
plane_weightings <- 6 / example_planes_df$quality

stereoplot(guides = FALSE)
points(planes, col = "grey", pch = 16, cex = .5)
contour(planes, add = TRUE, density.params = list(weights = plane_weightings))
```

![Diagram showing the density distribution some example data plotted in
an equal-area projection](Plots_files/figure-html/dens_contour-1.png)

[`contour()`](https://tobiste.github.io/structr/reference/stereo_contour.md)
adds contour lines, while
[`contourf()`](https://tobiste.github.io/structr/reference/stereo_contour.md)
adds filled contours and
[`image()`](https://tobiste.github.io/structr/reference/stereo_contour.md)
adds a density image (i.e. a dense grid of colored rectangles). See
[`?contour`](https://tobiste.github.io/structr/reference/stereo_contour.md),
[`?contourf`](https://tobiste.github.io/structr/reference/stereo_contour.md),
and
[`?image`](https://tobiste.github.io/structr/reference/stereo_contour.md)
for more information.

``` r
stereoplot(guides = FALSE)
points(planes, col = "grey", pch = 16, cex = .5)
contourf(planes, add = TRUE, density.params = list(weights = plane_weightings))
```

![Diagram showing the density distribution some example data plotted in
an equal-area projection](Plots_files/figure-html/dens_filled-1.png)

## Synopsis

Let’s create a publication ready synoptical plot for the line and plane
orientation data, showing the density distribution, eigenvectors/mean
values, and fabric strength.

``` r
# Minimum eigenvector of plane's orientation tensor:
planes_eigen <- ot_eigen(planes)$vectors
planes_eigen3 <- planes_eigen[3, ]

# Mean and SD of lines
lines_mean <- sph_mean(lines)
lines_sd <- sph_sd(lines)

# Fabric strength
fabric_p <- vollmer(planes)["D"]
fabric_l <- vollmer(lines)["D"]
```

The final plot:

``` r
# two plots side by side
par(mfrow = c(1, 2))

# Planes
stereoplot()
points(planes, col = "grey", pch = 16, cex = .5)
contour(planes, add = TRUE, weights = plane_weightings)
points(planes_eigen3, col = "black", pch = 16)
lines(planes_eigen3, col = "black", pch = 16)
title(
  main = "Planes",
  sub = paste0(
    "N: ", nrow(planes), " | Fabric strength: ", round(fabric_p, 2),
    "\nLambert equal area, lower hemisphere projection"
  )
)

# Lines
stereoplot()
points(lines, col = "grey", pch = 16, cex = .5)
contour(lines, add = TRUE, weights = line_weightings)
points(lines_mean, col = "black", pch = 16)
lines(lines_mean, ang = lines_sd, col = "black")
title(
  main = "Lines",
  sub = paste0(
    "N: ", nrow(lines), " | Fabric strength: ", round(fabric_l, 2),
    "\nLambert equal area, lower hemisphere projection"
  )
)
```

![Diagram showing the density distribution, eigenvectors and mean
values](Plots_files/figure-html/synops-1.png)

## Fault plots

Fault objects consist of planes (fault plane), lines (e.g. striae), and
the sense of movement. There are two ways how these combined features
can be visualized, namely the Angelier and the Hoeppener plot.

### Angelier plot

The **Angelier plot** shows all planes as *great circles* and lineations
as points (after Angelier, 1984)[⁶](#fn6). Fault striae are plotted as
vectors on top of the lineation pointing in the movement direction of
the hanging wall. Easy to read in case of homogeneous or small datasets.

``` r
f <- Fault(
  c("a" = 120, "b" = 125, "c" = 100),
  c(60, 62, 50),
  c(110, 25, 30),
  c(58, 9, 23),
  c(1, -1, 1)
)

stereoplot(title = "Angelier plot")
angelier(f, col = viridis::magma(nrow(f), end = .9))
```

![Diagram showing example fault data plotted in an equal-area
projection](Plots_files/figure-html/faults1-1.png)

### Hoeppener plot

The **Hoeppener plot** shows all planes as *poles* while lineations are
not shown (after Hoeppener, 1955)[⁷](#fn7). Instead, fault striae are
plotted as vectors on top of poles pointing in the movement direction of
the hanging wall. Useful in case of large or heterogeneous datasets.

``` r
stereoplot(title = "Hoeppener plot")
hoeppener(f, col = viridis::magma(nrow(f), end = .9), points = FALSE)
```

![Diagram showing example fault data plotted in an equal-area
projection](Plots_files/figure-html/faults2-1.png)

The `points` argument disables plotting the points at the start of the
arrows.

> [`fault_plot()`](https://tobiste.github.io/structr/reference/fault-plot.md)
> is a wrapper function that allows to switch between Angelier and
> Hoeppener plot using the `type` argument. See
> [`?fault_plot`](https://tobiste.github.io/structr/reference/fault-plot.md)
> for details.

## Assign plotting parameters based on data

{structr} offers some convenience functions to help your to map certain
plotting parameters such as color, size, and symbol based on vector
values:

- [`assign_col()`](https://tobiste.github.io/structr/reference/assign-color.md)
  for color

- [`assign_cex()`](https://tobiste.github.io/structr/reference/assign-cex.md)
  for marker size (character expansion)

- [`assign_pch()`](https://tobiste.github.io/structr/reference/assign-pch.md)
  for marker symbols (plotting character)

The following example assigns the `col` (color), the `cex` (size), and
the `pch` (symbol) based on three different properties:

``` r
# define three random properties
prop_continuous1 <- runif(nrow(planes))
prop_continuous2 <- runif(nrow(planes))
prop_discrete <- sort(letters[sample(1:3, size = nrow(planes), replace = TRUE)])

stereoplot()
points(planes, 
       col = assign_col(prop_continuous1), 
       cex = assign_cex(prop_continuous2, area = TRUE),
       pch = assign_pch(prop_discrete)
       )

# Add legends
legend_cex(prop_continuous2, area = TRUE, title = 'A continuous property', position = 'bottomleft', cex = .8)
legend_pch(prop_discrete, title = 'A discrete property', position = 'topleft', cex = .8)
legend_col(pretty(prop_continuous1), title = 'Another continuous property')
```

![](Plots_files/figure-html/assign-1.png)

The legend for these plotting parameters can be created using the
`legend_*` functions.

All these `assign_*` functions can be applied on continuous as well as
discrete values using `assign_*_d`. Also there are binned mapping
options through `assign_*_binnned`. See `?assign_cex()` for more
information.

------------------------------------------------------------------------

1.  Vollmer, F. W. (1990). An application of eigenvalue methods to
    structural domain analysis. Geological Society of America Bulletin,
    102, 786–791.

2.  Woodcock, N. H. (1977). Specification of fabric shapes using an
    eigenvalue method. Geological Society of America Bulletin88,
    1231–1236. Retrieved from
    <http://pubs.geoscienceworld.org/gsa/gsabulletin/article-pdf/88/9/1231/3418366/i0016-7606-88-9-1231.pdf>

3.  Lode, W. (1926). Versuche über den Einfluß der mittleren
    Hauptspannung auf das Fließen der Metalle Eisen. Kupfer und Nickel.
    Zeitschrift Für Physik, 36(11–12), 913–939.
    <https://doi.org/10.1007/BF01400222>

4.  Nádai, A. (1950). Theory of flow and fracture of solids.
    McGraw-Hill.

5.  Kamb, W. B. (1959). Ice Petrofabric Observations from Blue Glacier,
    Washington, in Relation to Theory and Experiment. Journal of
    Geophysical Research, 54(11).

6.  Angelier, J. Tectonic analysis of fault slip data sets, J. Geophys.
    Res. 89 (B7), 5835-5848 (1984)

7.  Hoeppener, R. Tektonik im Schiefergebirge. Geol Rundsch 44, 26-58
    (1955). <https://doi.org/10.1007/BF01802903>
