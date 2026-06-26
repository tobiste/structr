# Misfit parameters of slip inversion

Misfit parameters of slip inversion

## Usage

``` r
slip_inversion_misfit(tau, fault)
```

## Arguments

- tau:

  symmetric 3x3 matrix. The (reduced) stress tensor.

- fault:

  `"Fault"` object where the rows are the observations, and the columns
  the coordinates.

## Value

list.

- `"alpha"`:

  numeric. Angelier (1990)'s angles between the shear stress vector and
  the slip vector, ranging from 0&degree; (perfect fit) to 180°
  (inconsistent fit).

- `"beta"`:

  numeric. Michael (1984)'s angles between the tangential traction
  predicted by the best stress tensor and the slip vector on each plane,
  ranging from 0 to 90°.

- `"theta"`:

  numeric. Angle between slip planes and \\\sigma_1\\ ranging from 0 to
  180°.

- `"rup"`:

  numeric. "Ratio Upsilon" (RUP) parameter after Angelier (1990),
  ranging frm 0 (perfect fit) to 200% (misfit). See
  [`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md).

- `"quality"`:

  factor. Ranked misfit classification based on RUP values. See
  [`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md).

- `"quality_summary"`:

  integer. Counts of faults in the RUP-based quality ranks.

- `"misfit_means"`:

  Mean values of `alpha`, `beta`, `theta`, and `rup.`

## Examples

``` r
f <- angelier1990$TYM
tau <- reduced_stress(f)
slip_inversion_misfit(tau, f) 
#> $alpha
#>  [1] 19.0581833  0.4300616 17.6971171  8.5078919  9.3257870  5.9565593
#>  [7] 15.0070188  1.1696222 12.2657156  6.2970707 12.1518725  7.3901878
#> [13]  4.3095800 12.5543133  1.0365795 13.9582802  3.2573252  5.5898020
#> [19] 28.9782473 29.3190948 12.5920654 16.1272918 37.6709469  4.6992775
#> [25]  5.0514715 25.5265174 21.6148031 23.1509507  3.9315776 20.8077306
#> [31] 17.5255311  3.9282757 15.8194428 12.8556026  5.1748022  6.1157113
#> [37]  7.4467402 11.3905892
#> 
#> $beta
#>  [1] 22.4236622  7.9223563 20.1821988  9.9101003 13.9637462 20.1733189
#>  [7] 16.5519986  2.0706525 17.1911259  2.8583673 11.8477338  7.7629023
#> [13]  5.2680686 16.9325330  3.0864496 23.2285082  4.8017604 11.1360382
#> [19] 29.9292121 32.2496608 25.3837096 21.3747640 35.2452209 23.0365967
#> [25]  7.5504871 28.7097966 21.5712363 24.8114656  4.1798486  8.6411863
#> [31] 22.4333122  8.0510983 24.7064543 14.9681085  5.9747586  4.0508392
#> [37]  0.6771364 15.9340377
#> 
#> $theta
#>  [1] 61.29090 58.39885 75.84995 67.28394 63.01864 82.71051 75.13496 60.39790
#>  [9] 57.94764 58.44370 70.63671 51.01134 52.64142 65.43599 60.92119 62.12882
#> [17] 59.19148 62.47282 74.45211 53.74784 50.46237 43.95529 79.38046 44.13129
#> [25] 57.16206 76.39212 62.47759 55.78091 68.82664 56.60872 67.73859 62.87543
#> [33] 66.06326 45.75706 66.72009 70.86226 65.60097 67.09459
#> 
#> $rup
#>  [1] 39.266859 22.364271 42.794333 14.800803 20.461318 64.718606 38.120249
#>  [8] 18.789513 33.810029 26.242500 23.850645 35.834642 34.723973 22.583185
#> [15] 17.643757 28.502268 21.267580 16.202854 51.018881 66.347714 36.709113
#> [22] 50.411469 68.209016 34.135993 25.437568 50.886127 42.192872 53.871795
#> [29]  9.421844 47.329481 30.113931 13.797249 27.532282 46.592528  9.067210
#> [36] 20.237287 13.242937 19.749635
#> 
#> $quality
#>  [1] good       good       good       good       good       acceptable
#>  [7] good       good       good       good       good       good      
#> [13] good       good       good       good       good       good      
#> [19] acceptable acceptable good       acceptable acceptable good      
#> [25] good       acceptable good       acceptable good       good      
#> [31] good       good       good       good       good       good      
#> [37] good       good      
#> Levels: good acceptable poor
#> 
#> $misfit_means
#>    alpha     beta    theta      rup 
#> 12.22468 15.16869 62.65926 32.58638 
#> 
#> $quality_summary
#>       n_good n_acceptable       n_poor 
#>           31            7            0 
#> 
```
