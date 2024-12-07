
<!-- README.md is generated from README.Rmd. Please edit that file -->

# structr

<!-- badges: start -->

[![R-CMD-check](https://github.com/tobiste/structR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tobiste/structR/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

Structural geology package for R, free and open-source. It provides
functions to

- analyze and visualize orientation data of structural geology.

- analyze stress (including visualization of the magnitudes of stress in
  the Mohr Circle).

- reconstruct the orientation of structures in oriented drill cores
  using the alpha, beta, and gamma angles.

## Installation

You can install the development version of `{structr}` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tobiste/structr")
```

## Example

``` r
library(structr)
#> 
#> Attaching package: 'structr'
#> The following object is masked from 'package:stats':
#> 
#>     Pair
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.3.3

data(example_planes)
planes <- Plane(example_planes$dipdir, example_planes$dip)

fabric_strength <-  fabric_indexes(planes)

ggstereo() +
 geom_contourf_stereo(gg(planes), show.legend = TRUE) +
  geom_point(data = gg(planes), aes(x, y), color = 'lightgrey', shape = 21) +
  scale_fill_viridis_d(option = 'A' ) +
  labs(subtitle= "Example data", 
       title= 'Density contour plot',
       fill = 'Multiples of\nvon Mises-Fisher\ndistribution',
       caption = paste0('Equal area, lower hemisphere projection | Fabric strength: ', round(fabric_strength['C'], 2)))
#> Warning: Contour data has duplicated x, y coordinates.
#> ℹ 4900 duplicated rows have been dropped.
#> Warning: Removed 100 rows containing non-finite outside the scale range
#> (`stat_contour_filled()`).
```

<img src="man/figures/README-stereo-1.png" width="100%" />

## Documentation

in prep.

## Author

Tobias Stephan

## License

MIT License
