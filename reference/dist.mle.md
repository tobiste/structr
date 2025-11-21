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
#>        azimuth   plunge
#> mean  111.7604 60.96718
#> major 331.6261 23.07566
#> minor 234.3155 16.62996
#> 
#> $param
#>      kappa       beta        psi 
#>  3.6004056  0.4917274 -0.6043107 
#> 
#> $logcon
#> [1] 4.174607
#> 
#> $loglik
#> [1] -155.2405
#> 
vmf_mle(x)
#> $loglik
#> [1] -156.4088
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 111.76041  60.96718 
#> 
#> $kappa
#> [1] 3.550567
#> 
```
