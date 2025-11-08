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
#> mean  110.1078 60.48015
#> major 333.4227 22.39097
#> minor 235.6324 18.21170
#> 
#> $param
#>      kappa       beta        psi 
#>  3.6818192  0.5084907 -0.5752993 
#> 
#> $logcon
#> [1] 4.234482
#> 
#> $loglik
#> [1] -153.1613
#> 
vmf_mle(x)
#> $loglik
#> [1] -154.3944
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 110.10782  60.48015 
#> 
#> $kappa
#> [1] 3.625812
#> 
```
