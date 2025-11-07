# Oriented Drill Cores

This vignette describes how convert alpha-beta (and gamma) measurements
from oriented drill cores into geological planes and lines using the
[`drillcore_transformation()`](https://tobiste.github.io/structr/reference/drillcore.md)
function of the {structr} package.

``` r
library(structr)
```

Orientations in drill-cores are usually given by α and β angles
(lineations on a plane additionally have a γ angle) which describe
orientations with respect to the drill orientation. To convert these
angles from the “drillcore coordinate reference system” to our
geographical reference system, you can use the function
[`drillcore_transformation()`](https://tobiste.github.io/structr/reference/drillcore.md).

This function calculates the orientation of a plane or line from
internal core angles (**α**, **β**, and **γ**) of oriented drill cores,
using the the azimuth of drill core axis orientation (`azi` given in
degrees, measured clockwise from North), and the inclination of drill
core axis (`inc` in degrees, measured anticlockwise from the horizontal
plane).

``` r
azi <- 225
inc <- 45
```

> Note that negative values for the inclination indicate downward
> direction

`alpha` and `beta` are the internal core angles alpha and beta,
respectively, measured in degrees.

``` r
drillcore_transformation(azi, inc, alpha = 60, beta = 320)
#> Plane object (n = 1):
#> dip_direction           dip 
#>      25.00392      70.02959
```

The function returns a spherical objects. Since only alpha and beta
angles are specified, the output is a `"plane"` object.

For several alpha and beta angles:

``` r
planes_AB <- drillcore_transformation(azi, inc, alpha = c(60, 45), beta = c(320, 220))
```

The orientations can be plotted in a equal-area (lower hemisphere)
projection:

``` r
# initialize plot:
stereoplot()

# plot the core axis' azimuth and inclination)
points(Line(azi, inc))
text(Line(azi, inc), lab = "core-axis", pos = 1, font = 3)

# plot the plane orientations as poles...
points(planes_AB, col = c("#B63679FF", "#FB8861FF"))
text(planes_AB, lab = c("A", "B"), col = c("#B63679FF", "#FB8861FF"), pos = 1, font = 3)

# ... and as great circles
lines(planes_AB, col = c("#B63679FF", "#FB8861FF"))
```

![Diagram showing the orientation of the drillcore in a equal-area
projection](Oriented_Drill_Cores_files/figure-html/stereonet-1.png)

## References

Stigsson, M., & Munier, R. (2013). Orientation uncertainty goes bananas:
An algorithm to visualise the uncertainty sample space on stereonets for
oriented objects measured in boreholes. Computers and Geosciences, 56,
56–61. <https://doi.org/10.1016/j.cageo.2013.03.001>
