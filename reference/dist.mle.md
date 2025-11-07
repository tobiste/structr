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
#>        azimuth    plunge
#> mean  120.0725 53.289068
#> major 218.5818  6.296425
#> minor 313.1780 35.989104
#> 
#> $param
#>     kappa      beta       psi 
#> 4.3565562 0.6937126 0.6191050 
#> 
#> $logcon
#> [1] 4.749947
#> 
#> $loglik
#> [1] -137.3146
#> 
vmf_mle(x)
#> $loglik
#> [1] -139.3717
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 120.07246  53.28907 
#> 
#> $kappa
#> [1] 4.229803
#> 
```
