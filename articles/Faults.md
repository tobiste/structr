# Faults and Paleostress

This tutorial demonstrates the `"Fault"` object. It also shows how you
can extract paleo-stress directions from fault slip data, and how you
can derive displacement components from fault slip.

## Fault objects

A fault is given by the orientation of its plane (dip direction and dip
angle), the orientation of the slip (e.g. measured from striae, given in
azimuth and plunge angles), and the sense of displacement:

``` r

my_fault <- Fault(120, 50, 60, 110, sense = -1)
```

> Sense of fault displacement is 1 or -1 for normal or thrust offset,
> respectively.

Rake of the fault, i.e. the angle between fault slip vector and fault
strike:

``` r

Fault_rake(my_fault)
```

    ## [1] -107.2294

Define a fault by just knowing the orientation of the fault plane, the
sense, and the rake

``` r

# 1. Define a plane through dip direction, dip angle
fault_plane <- Plane(c(120, 120, 100), c(60, 60, 50))

# 2. Define a fault through the plane and rake angle:
Fault_from_rake(fault_plane, rake = c(84.7202, -10, 30))
```

    ## Fault object (n = 3):
    ##      dip_direction dip   azimuth    plunge sense
    ## [1,]           120  60 109.52858 59.581591     1
    ## [2,]           120  60  24.96163 -8.649165    -1
    ## [3,]           100  50  30.36057 22.521012     1

Often, measured orientation angles can be (slightly) imprecise and
subjected to some random noise. Thus the slip vector will not lie
(perfectly) on the fault plane, judging by the measurements. To correct
the measurements so that this will not be the case:

``` r

p <- Pair(120, 60, 110, 58)
```

``` r

misfit_pair(p)
```

    ## $fvec
    ## Vector (Vec3) object (n = 1):
    ##          x          y          z 
    ##  0.4306083 -0.7432653  0.5119893 
    ## 
    ## $lvec
    ## Vector (Vec3) object (n = 1):
    ##          x          y          z 
    ## -0.1752490  0.4876331  0.8552787 
    ## 
    ## $misfit
    ## [1] 0.02793105

``` r

correct_pair(p)
```

    ## Pair object (n = 1):
    ## dip_direction           dip       azimuth        plunge 
    ##     120.08572      59.20357     109.76778      58.79054

> A `"Pair"` object is a container of associated plane and line
> measurements. Basically like a fault without the sense of
> displacement.

## Plotting faults

Fault objects consist of planes (fault plane), lines (e.g. striae), and
the sense of movement. There are two ways how these combined features
can be visualized, namely the Angelier and the Hoeppener plot.

- The **Angelier plot** shows all planes as *great circles* and
  lineations as points (after Angelier, 1984). Fault striae are plotted
  as vectors on top of the lineation pointing in the movement direction
  of the hanging wall. Easy to read in case of homogeneous or small
  datasets.
- The **Hoeppener plot** shows all planes as *poles* while lineations
  are not shown (after Hoeppener, 1955). Instead, fault striae are
  plotted as vectors on top of poles pointing in the movement direction
  of the hanging wall. Useful in case of large or heterogeneous
  datasets.

``` r

# simongomez is a example fault dataset:

# define some colors for each fault in the dataset (here the fault sense)
fault_cols <- assign_col(simongomez[, 5], pal = viridis::magma, begin = .2, end = .8)

par(mfrow = c(1, 2))
stereoplot(title = "Angelier plot")
angelier(simongomez, col = fault_cols)

stereoplot(title = "Hoeppener plot")
hoeppener(simongomez, col = fault_cols, points = FALSE)
```

![Diagram showing some fault data in equa-area projection using the
Angelier and Hoeppener representation of
fault-slip](Faults_files/figure-html/faults_plotting-1.png)

## Fault stress analysis

The *Wallace-Bott hypothesis* states that fault slip occurs parallel to
the maximum shear stress. This allows to reconstruct stress axes using
fault-slip data. {structr} offers several techniques to calculate the
orientation of principal stress axes, the simple P-T method, and a
fault-slip inversion technique.

### P-T method

This simple technique calculates PT-axes, kinematic plane (M), and
dihedra separation plane (d).

First we load some example data (here the first three faults from the
TYM dataset by from Angelier, 1990)[^1]

``` r

data("angelier1990")
my_fault2 <- angelier1990$TYM[1:3, ]

print(my_fault2)
```

    ## Fault object (n = 3):
    ##      dip_direction dip  azimuth   plunge sense
    ## [1,]           137  61 117.0135 59.46650     1
    ## [2,]           128  59 146.8990 57.58045     1
    ## [3,]             2  80 287.5304 56.63295     1

``` r

my_fault2_PT <- Fault_PT(my_fault2)
print(my_fault2_PT)
```

    ## $p
    ## Line object (n = 3):
    ##       azimuth   plunge
    ## [1,] 340.6153 72.15071
    ## [2,] 281.6092 73.96385
    ## [3,] 214.3217 45.50717
    ## 
    ## $t
    ## Line object (n = 3):
    ##       azimuth   plunge
    ## [1,] 129.6815 15.44075
    ## [2,] 135.2532 13.45689
    ## [3,] 336.9158 27.88916
    ## 
    ## $m
    ## Plane object (n = 3):
    ##      dip_direction      dip
    ## [1,]      42.11395 81.26433
    ## [2,]     223.18909 81.43997
    ## [3,]     265.80721 58.54232
    ## 
    ## $d
    ## Plane object (n = 3):
    ##      dip_direction      dip
    ## [1,]      297.0135 30.53350
    ## [2,]      326.8990 32.41955
    ## [3,]      107.5304 33.36705

Plot the results

``` r

stereoplot(title = "PT results", guides = FALSE)
fault_plot(my_fault2)
points(my_fault2_PT$p, col = "#B63679FF", pch = 16)
points(my_fault2_PT$t, col = "#FEC287FF", pch = 18)
lines(my_fault2_PT$t, lty = 2, col = "grey40")
lines(my_fault2_PT$d, lty = 3, col = "grey80")

legend("right",
  legend = c("P-axis", "T-axis", "M-plane", "Diheadra"),
  col = c("#B63679FF", "#FEC287FF", "grey40", "grey80"),
  pch = c(16, 18, NA, NA), lty = c(NA, NA, 2, 3)
)
```

![Diagram showing some fault-slip data and the PT axes in a
stereoplot](Faults_files/figure-html/faults_pt_plot-1.png)

### Fault slip inversion

Our goal is to find the single uniform stress tensor that most likely
caused the faulting events. With only slip data to constrain the stress
tensor the isotropic component can not be determined, unless assumptions
about the fracture criterion are made. Hence inversion will be for the
deviatoric stress tensor only. A single fault can not completely
constrain the deviatoric stress tensor a, therefore it is necessary to
simultaneously solve for a number of faults, so that a single a that
best satisfies all of the faults is found.

> This is equivalent to assuming that the stress field is a constant
> tensor within the region being studied for the duration of the
> faulting event.

{structr} provides two numerical solutions to determine the orientation
of the principal stresses from fault slip data.

- Michael (1984) linear inversion method[^2] which uses bootstrapping
  for confidence intervals of the stress estimates.

- Angelier (1990) direct inversion method[^3] coupled with the iterative
  optimization after Mostafa (2005)[^4] to find the best fit reduced
  stress tensor.

First we load some example data (here the data from Angelier, 1990)[^5]

``` r

test_data <- angelier1990$TYM

stereoplot(title = "Test data", guides = FALSE)
fault_plot(test_data, col = "grey30")
```

![Diagram showing some fault-slip data in a
stereoplot](Faults_files/figure-html/slip_inversion1-1.png)

The stress inversion using the Michael method with 10 bootstraps:

``` r

test_res <- slip_inversion(test_data, method = "michael", n_iter = 10)

# Average beta angle
test_res$misfit$beta
```

    ##  [1] 22.4236622  7.9223563 20.1821988  9.9101003 13.9637462 20.1733189
    ##  [7] 16.5519986  2.0706525 17.1911259  2.8583673 11.8477338  7.7629023
    ## [13]  5.2680686 16.9325330  3.0864496 23.2285082  4.8017604 11.1360382
    ## [19] 29.9292121 32.2496608 25.3837096 21.3747640 35.2452209 23.0365967
    ## [25]  7.5504871 28.7097966 21.5712363 24.8114656  4.1798486  8.6411863
    ## [31] 22.4333122  8.0510983 24.7064543 14.9681085  5.9747586  4.0508392
    ## [37]  0.6771364 15.9340377

``` r

# Average resolved shear stress
test_res$tau_mean
```

    ## [1] 0.9289024

To check the accuracy of the solution, you can evaluate

- the average angle β between the tangential traction predicted by the
  best stress tensor and the slip vector on each plane (ideally close to
  0), and
- the average resolved shear stress on each plane (should be close to
  1).

To visualizing the orientation of the principal stresses and the
confidence region of the axes, you may use the stereoplot() functions
from {structr}

``` r

cols <- c("#000004FF", "#B63679FF", "#FEC287FF")

stereoplot(title = "Stress inversion", guides = FALSE)
fault_plot(test_data, col = "grey75")
stereo_confidence(test_res$principal_axes_CI$sigma1, col = cols[1])
stereo_confidence(test_res$principal_axes_CI$sigma2, col = cols[2])
stereo_confidence(test_res$principal_axes_CI$sigma3, col = cols[3])
text(test_res$principal_axes,
  label = rownames(test_res$principal_axes),
  col = cols, adj = -.25
)
legend("topleft",
  col = cols,
  legend = rownames(test_res$principal_axes), pch = 16
)
```

![Diagram showing principal stress vector results in a
stereoplot](Faults_files/figure-html/slip_inversion_plot-1.png)

The stress shape ratio Φ (Angelier 1979)[^6]

``` r

test_res$stress_shape$phi
```

    ##   sigma2 
    ## 0.101247

``` r

# 95% confidence interval
test_res$phi_CI
```

    ## [1] 0.08715034 0.15052106
    ## attr(,"conf.level")
    ## [1] 0.95

The angle β is the angle between the tangential traction predicted by
the best stress tensor and the slip vector. This deviation can be
visualized in the stereoplot:

``` r

beta <- test_res$misfit$beta

stereoplot(
  title = "Deviation",
  sub = bquote(bar(beta) == .(round(test_res$beta)) * degree),
  guides = FALSE
)

fault_plot(test_data, col = assign_col(beta))
legend_col(
  seq(min(beta), max(beta), 10),
  title = bquote("Deviation angle" ~ beta ~ "(" * degree * ")")
)
```

![Diagram showing deviation of the best-fit principal stress tensor in a
stereoplot](Faults_files/figure-html/slip_inversion_beta-1.png)

``` r

Mohr_plot(
  sigma1 = test_res$principal_vals[1],
  sigma2 = test_res$principal_vals[2],
  sigma3 = test_res$principal_vals[3],
  unit = NULL, include.zero = FALSE
)
points(test_res$stress_components[, 'normal'], abs(test_res$stress_components[, 'shear']),
  col = assign_col(beta), pch = 16
)
```

![Diagram showing the Mohr circle, and the normal and shear stresses for
the slip inversion
result](Faults_files/figure-html/slip_inversion_mohr-1.png)

### Maximum horizontal stress

The orientation of the maximum horizontal stress
($`\sigma_\text{Hmax}`$) can be calculated from the stress tensor the
the orientation of the principal stress ($`\sigma_1`$, $`\sigma_2`$,
$`\sigma_3`$) axes their their relative magnitudes($`R`$) [^7].

First, we define the orientation of the principle stress axes:

``` r

S1 <- Line(250.89, 70.07)
S3 <- Line(103.01, 17.07)
```

To get $`\sigma_2`$ - which is perpendicular to $`\sigma_1`$ and
$`\sigma_3`$, we calculate the cross-product of the two vectors:

``` r

S2 <- crossprod(S3, S1)
```

The azimuth of $`\sigma_\text{Hmax}`$ for a given stress ratio `R = 1`:

``` r

SH(S1, S2, S3, R = 1) # in degrees
```

    ## [1] 70.89

For a several stress ratios:

``` r

R <- seq(0, 1, .1)
cbind(R, SH = SH(S1, S2, S3, R = R))
```

    ##         R       SH
    ##  [1,] 0.0 13.01021
    ##  [2,] 0.1 13.37695
    ##  [3,] 0.2 13.84162
    ##  [4,] 0.3 14.44908
    ##  [5,] 0.4 15.27621
    ##  [6,] 0.5 16.46586
    ##  [7,] 0.6 18.31445
    ##  [8,] 0.7 21.53704
    ##  [9,] 0.8 28.23884
    ## [10,] 0.9 45.01043
    ## [11,] 1.0 70.89000

The $`\sigma_\text{Hmax}`$ for our slip inversion result from above:

``` r

SH(
  S1 = test_res$principal_axes[1, ],
  S2 = test_res$principal_axes[2, ],
  S3 = test_res$principal_axes[3, ],
  R = test_res$stress_shape$R
)
```

    ## [1] 60.80844

## Fault offsets

The offset along a fault can be factorized into several components.

![Fig. 1: Graphic illustration of displacement components along a fault.
\mathbb{G} as the georeference frame with D = down, E = East, N =
North.](fault_displacements.png)

Fig. 1: Graphic illustration of displacement components along a fault.
$`\mathbb{G}`$ as the georeference frame with D = down, E = East, N =
North.

### Get different components with trigonometry

#### 1. Input: Fault orientation (dip angle, dip direction), shortening direction, and horizontal throw

Knowing the horizontal throw (e.g. from plate motion parameters), the
remaining components of the displacements along a given fault are as
follows.

``` math
\begin{equation}\delta = |\sigma_{\textrm{Hmax}} - (\textrm{dip direction}+90^{\circ})|\end{equation}
```

*Slip components in the horizontal plane:*

``` math
\begin{equation}\begin{split}
  f_\textrm{strike slip} & = |\cos{\delta} * f_\textrm{horizontal throw}|\\
  f_\textrm{heave} & = \sqrt{f_\textrm{horizontal throw}^2 - f_\textrm{strike slip}^2}
  \end{split}\end{equation}
```

*Slip components in the vertical plane perpendicular to the strike of
the fault:*

``` math
\begin{equation}\begin{split}
  f_\textrm{dip slip} &=  \frac{f_\textrm{heave}}{cos{(\textrm{dip})}}\\
  f_\textrm{vertical throw} &= \sqrt{f_\textrm{dip slip}^2 - f_\textrm{heave}^2}
  \end{split}\end{equation}
```

*Slip components in the fault plane plane:*

``` math
\begin{equation}\begin{split}
  f_\textrm{net slip} &= \sqrt{f_\textrm{strike slip}^2 + f_\textrm{dip slip}^2}\\
  \textrm{rake} &= \arctan{\left(\frac{f_\textrm{dip slip}}{f_\textrm{strike slip}}\right)}
  \end{split}\end{equation}
```

Thus, the rake angle describes the ratio between the dip slip and the
strike slip component.

Knowing the vertical throw (e.g. from thermochronology or petrology),
the fault dip (?assumption), and the direction and amount of horizontal
offset, the **strike** of the fault is as follows:

``` math
\begin{equation}\begin{split}
f_\textrm{heave} &= f_\textrm{vertical throw} * \tan{(\textrm{dip})}\\
\delta &= \arcsin{\left(\frac{f_\textrm{heave}}{f_\textrm{horizontal throw}}\right)}\\
\textrm{strike} &= |\sigma_{\textrm{Hmax}} - \delta|
\end{split}\end{equation}
```

Knowing the vertical throw (e.g. from thermochronology or petrology),
the fault strike (geomorphology), and the direction and amount of
horizontal offset, the **dip** of the fault is as follows:

``` math
\begin{equation}\begin{split}
\delta &= |\sigma_{\textrm{Hmax}} - \textrm{strike}|\\
f_\textrm{heave} &= f_\textrm{horizontal throw} * \sin{\delta} \\
\textrm{dip} &= \arctan{\left(\frac{f_\textrm{heave}}{f_\textrm{vertical throw}}\right)} 
\end{split}\end{equation}
```

Knowing the vertical throw (e.g. from thermochronology or petrology) and
the fault’s dip and rake, the *horizontal offset*, the horizontal throw,
and the net-slip are as follows:

### Fault displacement tensors

Each fault component is a vector describing its direction and length.
For instance, the vector of the strike slip is:

``` math
\begin{equation}
  \vec{f_\text{strike slip}} = \begin{pmatrix} \lVert f_\textrm{strike-slip}\rVert \\ 0 \\ 0 \end{pmatrix}
\end{equation}
```

### Net slip vector

Net slip vector

``` math
\begin{equation}
  \vec{f_\text{net}} = \begin{pmatrix} \lVert f_\textrm{strike-slip}\rVert \\ \lVert f_\textrm{heave}\rVert \\ \lVert f_\textrm{vertical throw}\rVert \end{pmatrix}
\end{equation}
```

``` r

fault_displacements(strikeslip = 2, verticalthrow = -5, heave = 3)
```

### Principal displacement tensor

The Eigen values of $`\mathsf{F}_{ij}`$; represented as
$`\{ \vec{f_1}, \vec{f_2}, \vec{f_3} \}`$ are referred to as the heave,
strike slip, and vertical throw component, respectively. These
orthonormal vectors define a orthogonal matrix, i.e. the **principal
displacement tensor** $`\mathsf{F}_\mathbb{F}`$:

``` math
\begin{equation}\mathsf{F}_{\mathbb{F}} = {\begin{bmatrix} \lVert f_\textrm{heave}\rVert & 0 & 0\\ 0 & \lVert f_\textrm{strike-slip}\rVert & 0\\ 0 & 0 & \lVert f_\textrm{vertical throw}\rVert\end{bmatrix}}\end{equation}
```

The tensor $`\mathsf{F}`$ can also be defined by the magnitudes of the
fault displacements:

``` r

Fu <- displacement_tensor(s = 2, v = -5, h = 3)
print(Fu)
```

    ##      [,1] [,2] [,3]
    ## [1,]    3    0    0
    ## [2,]    0    2    0
    ## [3,]    0    0   -5
    ## attr(,"class")
    ## [1] "matrix"  "array"   "ftensor"

### Orientation tensor

Fault orientation tensor is defined by the fault plane’s location,
orientation (dip direction and dip angle), and the fault’s slip
(direction and magnitude):

``` math
\begin{equation}\mathsf{F}_{ij} = {\begin{bmatrix}f_{11} & f_{12} & f_{13}\\ f_{21} & f_{22} & f_{23}\\ f_{31} & f_{32} & f_{33}\end{bmatrix}}\end{equation}
```

``` r

Fg <- displacement_tensor(s = 2, v = -5, h = 3, dip_direction = 45)
print(Fg)
```

    ##         [,1]      [,2] [,3]
    ## [1,] 2.12132  1.414214    0
    ## [2,] 2.12132 -1.414214    0
    ## [3,] 0.00000  0.000000   -5
    ## attr(,"class")
    ## [1] "matrix"  "array"   "ftensor"

### From Principal displacement tensor to Orientation tensor

Translation point of origin in $`\mathsf{F_\mathbb{F}}`$ into point of
measurement and rotate into fault orientation $`\mathsf{F}_\mathbb{FG}`$

``` r

displacement_tensor_decomposition(Fg, dip_direction = 45)
```

    ## $displacements
    ##           dip    delta     rake verticalthrow horizontalthrow heave  dipslip
    ## [1,] 59.03624 56.30993 71.06818            -5        3.605551     3 5.830952
    ##      strikeslip  netslip
    ## [1,]          2 6.164414
    ## 
    ## $fault
    ## Fault object (n = 1):
    ## dip_direction           dip       azimuth        plunge         sense 
    ##      45.00000      59.03624     258.69007     -54.20424      -1.00000 
    ## 
    ## $strain_tensor
    ##      [,1] [,2] [,3]
    ## [1,]  3.0 -1.5 -1.5
    ## [2,] -1.0  1.0 -1.0
    ## [3,]  2.5  2.5 15.0
    ## 
    ## $volumetric_strain
    ## [1] 19
    ## 
    ## $shear_strain
    ## [1] 5.196152

## References

Angelier, J. (1979). Determination of the mean principal directions of
stresses for a given fault population. Tectonophysics, 56(3–4), T17–T26.
<https://doi.org/10.1016/0040-1951(79)90081-7>

Angelier, J. (1990). Inversion of field data in fault tectonics to
obtain the regional stress—III. A new rapid direct inversion method by
analytical means. Geophys. J. Int, 103, 363–376.
<https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

Lund, B., & Townend, J. (2007). Calculating horizontal stress
orientations with full or partial knowledge of the tectonic stress
tensor. Geophysical Journal International, 170(3), 1328–1335.
<https://doi.org/10.1111/j.1365-246X.2007.03468.x>

Michael, A. J. (1984). Determination of stress from slip data: Faults
and folds. Journal of Geophysical Research: Solid Earth, 89(B13),
11517–11526. <https://doi.org/10.1029/JB089iB13p11517>

Mostafa, M. E. (2005). Iterative direct inversion: An exact
complementary solution for inverting fault-slip data to obtain
palaeostresses. Computers & Geosciences, 31(8), 1059–1070.
<https://doi.org/10.1016/j.cageo.2005.02.012>

[^1]: Angelier, J. (1990). Inversion of field data in fault tectonics to
    obtain the regional stress—III. A new rapid direct inversion method
    by analytical means. Geophys. J. Int, 103, 363–376.
    <https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

[^2]: Michael, A. J. (1984). Determination of stress from slip data:
    Faults and folds. Journal of Geophysical Research: Solid Earth,
    89(B13), 11517–11526. <https://doi.org/10.1029/JB089iB13p11517>

[^3]: Angelier, J. (1990). Inversion of field data in fault tectonics to
    obtain the regional stress—III. A new rapid direct inversion method
    by analytical means. Geophys. J. Int, 103, 363–376.
    <https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

[^4]: Mostafa, M. E. (2005). Iterative direct inversion: An exact
    complementary solution for inverting fault-slip data to obtain
    palaeostresses. Computers & Geosciences, 31(8), 1059–1070.
    <https://doi.org/10.1016/j.cageo.2005.02.012>

[^5]: Angelier, J. (1990). Inversion of field data in fault tectonics to
    obtain the regional stress—III. A new rapid direct inversion method
    by analytical means. Geophys. J. Int, 103, 363–376.
    <https://doi.org/10.1111/j.1365-246X.1990.tb01777.x>

[^6]: Angelier, J. (1979). Determination of the mean principal
    directions of stresses for a given fault population. Tectonophysics,
    56(3–4), T17–T26. <https://doi.org/10.1016/0040-1951(79)90081-7>

[^7]: Lund & Townend (2007): Calculating horizontal stress orientations
    with full or partial knowledge of the tectonic stress tensor.
    *Geophys. J. Int.*, 170, 1328—1335. doi:
    10.1111/j.1365-246X.2007.03468.x
