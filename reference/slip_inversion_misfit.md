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

- `alpha`:

  numeric. Deviation angle between slickenline and shear traction.
  Ranging from 0&degree; (perfect fit) to 90° (inconsistent fit).

- `alpha_mean`:

  numeric. The mean of `alpha`, the ie. the mean deviation of predicted
  from observed slip.

- `rup`:

  numeric. "Ratio Upsilon" (RUP) parameter after Angelier (1990),
  ranging frm 0 (perfect fit) to 200% (misfit). See
  [`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md).

- `quality`:

  factor. Ranked misfit classification based on RUP values. See
  [`tau2rup()`](https://tobiste.github.io/structr/reference/tau2rup.md).

- `rup_mean`:

  numeric. The mean RUP.

- `quality_summary`:

  integer. Counts of faults in the RUP-based quality ranks.

- `flipped`:

  logical. Are the signs of slip vectors flipped, i.e., dot product of
  the slip ray and the predicted shear traction is negative?

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
#> $alpha_mean
#> [1] 12.25499
#> 
#> $rup
#>  [1] 39.266859 22.364271 42.794333 14.800803 20.461318 64.718606 38.120249
#>  [8] 18.789513 33.810029 26.242500 23.850645 35.834642 34.723973 22.583185
#> [15] 17.643757 28.502268 21.267580 16.202854 51.018881 66.347714 36.709113
#> [22] 50.411469 68.209016 34.135993 25.437568 50.886127 42.192872 53.871795
#> [29]  9.421844 47.329481 30.113931 13.797249 27.532282 46.592528  9.067210
#> [36] 20.237287 13.242937 19.749635
#> 
#> $rup_mean
#> [1] 32.58638
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
#> $quality_summary
#>       n_good n_acceptable       n_poor 
#>           31            7            0 
#> 
#> $flipped
#>  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [25] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [37] FALSE FALSE
#> 
```
