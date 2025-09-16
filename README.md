
<!-- README.md is generated from README.Rmd. Please edit that file -->

# structr

<!-- badges: start -->

[![R-CMD-check](https://github.com/tobiste/structr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tobiste/structr/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`{structr}` is a free and open-source package for R that provides tools
for structural geology. You can

- analyze and visualize orientation data of structural geology
  (including, stereographic projections, contouring, fabric plots, and
  statistics),

- do a cluster analysis of orientation data,

- analyze stress (including visualization of the magnitudes of stress in
  the **Mohr circle** and extracting the maximum horizontal stress of a
  3D stress tensor),

- reconstruct the orientation of structures in oriented drill cores
  using the alpha, beta, and gamma angles,

- calculate fault displacement components,

- reconstruct orientation and magnitudes using fault-slip inversion
  (**Paleostress analysis**),

- contour geologic fabric and finite strain data (**R**<sub>f</sub>/ϕ)
  on the unit hyperboloid,

- perform vorticity analysis using the **Rigid Grain Net** method, and

- directly import your field data from **StraboSpot**.

> The {structr} package is all about structures in 3D. For analyzing
> orientations in 2D (statistics, rose diagrams, etc.), check out the
> [tectonicr](https://github.com/tobiste/tectonicr) package!

## Installation

You can install the development version of `{structr}` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tobiste/structr")
```

## Example

This is a basic example which shows you how to import data and plot them
on stereographic or equal-area projections.

``` r
library(structr)
library(ggplot2)
library(mapproj)

data(example_planes)
planes <- Plane(example_planes$dipdir, example_planes$dip)

fabric <- or_shape_params(planes)$Vollmer["D"]

ggstereo() +
  geom_contourf_stereo(gg(planes)) +
  geom_point(data = gg(planes), aes(x, y), color = "grey", shape = 21) +
  scale_fill_viridis_d(option = "F") +
  labs(
    subtitle = "Example data",
    title = "Density contour plot",
    fill = "Multiples of\nvon Mises-Fisher\ndistribution",
    caption = paste0("Equal area, lower hemisphere projection | N: ", nrow(planes), " | Fabric strength: ", round(fabric, 2))
  )
```

<img src="man/figures/README-stereo-1.png" width="100%" />

## Documentation

The detailed documentation can be found at
<https://tobiste.github.io/structr/>

## Author

Tobias Stephan (<tstephan@lakeheadu.ca>)

## Feedback, issues, and contributions

I welcome feedback, suggestions, issues, and contributions! If you have
found a bug, please file it
[here](https://github.com/tobiste/structr/issues) with minimal code to
reproduce the issue.

## License

MIT License
