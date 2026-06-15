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
#>         azimuth   plunge
#> mean  124.68842 54.04522
#> major  24.39746  7.38321
#> minor 289.20073 34.95350
#> 
#> $param
#>     kappa      beta       psi 
#> 3.8607553 0.5336659 0.2900435 
#> 
#> $logcon
#> [1] 4.366946
#> 
#> $loglik
#> [1] -148.6928
#> 
vmf_mle(x)
#> $loglik
#> [1] -149.9883
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 124.68842  54.04522 
#> 
#> $kappa
#> [1] 3.794865
#> 
```
