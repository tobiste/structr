# MLE of spherical rotational symmetric distributions

Estimates the parameters of a von Mises-Fisher or Kent distribution.

## Usage

``` r
kent_mle(x)

vmf_mle(x)
```

## Source

Adapted from
[`Directional::kent.mle()`](https://rdrr.io/pkg/Directional/man/kent.mle.html)
and `Directional::vmf.mle()`

## Arguments

- x:

  object of class `"Vec3"`, `"Line"`, `"Ray"`, `"Plane"`, `"Pair"`, or
  `"Fault"`, where the rows are the observations and the columns are the
  coordinates.

## Examples

``` r
x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
kent_mle(x)
#> $G
#> Line object (n = 3):
#>         azimuth    plunge
#> mean  119.66984 57.805223
#> major  24.20474  3.431543
#> minor 292.06018 31.966699
#> 
#> $param
#>     kappa      beta       psi 
#> 3.7419858 0.5657010 0.3366434 
#> 
#> $logcon
#> [1] 4.281925
#> 
#> $loglik
#> [1] -151.6588
#> 
vmf_mle(x)
#> $loglik
#> [1] -153.2303
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 119.66984  57.80522 
#> 
#> $kappa
#> [1] 3.66987
#> 
```
