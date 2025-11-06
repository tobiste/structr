# Strain and Vorticity Analysis

This tutorial demonstrates how to perform strain and vorticity analysis
using the {structr} package in R.

``` r
library(structr)
```

## Strain analysis

Import some Rf/ϕ data from elliptical strain markers

``` r
data(ramsay)
head(ramsay)
#>         R   phi
#> [1,] 1.24 35.96
#> [2,] 1.52 27.59
#> [3,] 1.33 36.91
#> [4,] 1.78 27.31
#> [5,] 1.51 17.73
#> [6,] 1.70 24.45
```

The mean mean strain ellipse (shape and orientation) of deformed
elliptical objects as strain markers can be calculated by using the mean
shape matrix and its eigenvalues[¹](#fn1):

``` r
ramsay_mean <- mean_strain_ellipse(r = ramsay[, 1], phi = ramsay[, 2])
print(ramsay_mean)
#> $R
#> [1] 1.628138
#> 
#> $phi
#> [1] 25.73632
#> 
#> $R_CI
#> [1] 1.59275 1.66359
#> 
#> $phi_CI
#> [1] 24.71905 26.75180
```

The {structr} algorithm also calculates bootstrapped 95% confidence
interval.

To visualize the distribution of the strain values, we can calculate
densities in Rf/ϕ space[²](#fn2) and plot them in a **Rf/ϕ**
diagram[³](#fn3):

``` r
Rphi_plot(r = ramsay[, 1], phi = ramsay[, 2])
```

![Rf-phi plot](Strain_Analysis_files/figure-html/rphi-1.png)

or in an **Equidistant polar plot**[⁴](#fn4):

``` r
Rphi_polar_plot(ramsay[, 1], ramsay[, 2], proj = "eqd")
```

![R/phi data in polar
plot](Strain_Analysis_files/figure-html/plot2-1.png)

## 3D Strain

Three-dimensional strain data are expressed by the ratio of the
magnitudes of the 3 principal strain axes of the strain ellipsoid,
$X \geq Y \geq Z$. They can be represented in the **Flinn
diagram**[⁵](#fn5), either linear or logarithmic.

``` r
data("holst")
R_XY <- holst[, "R_XY"]
R_YZ <- holst[, "R_YZ"]

flinn_plot(cbind(R_XY, R_YZ), log = TRUE, col = "#B63679", pch = 16)
```

![Flinn diagram](Strain_Analysis_files/figure-html/flinn-1.png)

or the **Hsu diagram**[⁶](#fn6) using the natural octahedral unit
strain[⁷](#fn7) (${\bar{\epsilon}}_{s}$) and the Lode parameter[⁸](#fn8)
($\nu$):

``` r
hsu_plot(cbind(R_XY, R_YZ), col = "#B63679", pch = 16)
```

![Hsu diagram](Strain_Analysis_files/figure-html/hsu-1.png)

## Vorticity analysis

The **rigid grain net** after (Jessup et al. 2007)[⁹](#fn9) plots the
distribution the strain ratio (`R`) of orientation (`phi`) of
porphyroclast over the theoretical distribution of tailless clasts. The
plot estimates the critical threshold `Rc` marking the transition
between the stable-sink position and infinitely rotating porphyroclasts.
This threshold can be interpreted as the the **mean kinematic vorticity
number**. Here the `Rc` is estimated using the bootstrap method
described in Stephan et al. (2025)[¹⁰](#fn10).

``` r
data(shebandowan)
set.seed(20250411)

# Color code porphyroclasts by size of clast (area in log-scale):
RGN_plot(shebandowan$r, shebandowan$phi, col = assign_col(log(shebandowan$area)), pch = 16)
```

![Rigid grain net](Strain_Analysis_files/figure-html/vorticity-1.png)

------------------------------------------------------------------------

1.  Shimamoto, T., & Ikeda, Y. (1976). A simple algebraic method for
    strain estimation from deformed ellipsoidal objects. 1. Basic
    theory. Tectonophysics, 36(4), 315–337.
    <https://doi.org/10.1016/0040-1951(76)90107-4>

2.  Vollmer, F. W. (2018). Computers and Geosciences Automatic
    contouring of geologic fabric and finite strain data on the unit
    hyperboloid. Computers and Geosciences, 115(June 2017), 134–142.
    <https://doi.org/10.1016/j.cageo.2018.03.006>

3.  Ramsay, J. G., & Huber, M. I. (1983). The Techniques of Modern
    Structural Geology: Strain Analyses (Vol. 1). London: Academic
    Press.

4.  Elliott, D. (1970). Determination of Finite Strain and Initial Shape
    from Deformed Elliptical Objects. GSA Bulletin, 81(8), 2221–2236.
    <https://doi.org/https://doi.org/10.1130/0016-7606(1970)81%5B2221:DOFSAI%5D2.0.CO;2>

5.  Flinn, D. (1965). On the Symmetry Principle and the Deformation
    Ellipsoid. Geological Magazine, 102(1), 36–45.
    <https://doi.org/10.1017/S0016756800053851>

6.  Hsu, T. C. (1966). The characteristics of coaxial and non-coaxial
    strain paths. Journal of Strain Analysis, 1(3), 216–222.
    <https://doi.org/10.1243/03093247V013216>

7.  Nádai, A. (1950). Theory of flow and fracture of solids.
    McGraw-Hill.

8.  Lode, W. (1926). Versuche über den Einfluß der mittleren
    Hauptspannung auf das Fließen der Metalle Eisen. Kupfer und Nickel.
    Zeitschrift Für Physik, 36(11–12), 913–939.
    <https://doi.org/10.1007/BF01400222>

9.  Jessup, M. J., Law, R. D., & Frassi, C. (2007). The Rigid Grain Net
    (RGN): An alternative method for estimating mean kinematic vorticity
    number (Wm). Journal of Structural Geology, 29(3), 411–421.
    <https://doi.org/10.1016/j.jsg.2006.11.003>

10. Stephan, T., Phillips, N., Tiitto, H., Perez, A., Nwakanma, M.,
    Creaser, R., & Hollings, P. (2025). Going with the flow — Changes of
    vorticity control gold enrichment in Archean shear zones
    (Shebandowan Greenstone Belt, Superior Province, Canada). Journal of
    Structural Geology, 201(September), 105542.
    <https://doi.org/10.1016/j.jsg.2025.105542>
