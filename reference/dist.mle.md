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
#> mean  118.6167 56.19982
#> major 358.3694 18.63526
#> minor 258.4373 27.08863
#> 
#> $param
#>      kappa       beta        psi 
#>  4.3791988  0.2902575 -0.1862201 
#> 
#> $logcon
#> [1] 4.747578
#> 
#> $loglik
#> [1] -136.4245
#> 
vmf_mle(x)
#> $loglik
#> [1] -136.4999
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 118.61668  56.19982 
#> 
#> $kappa
#> [1] 4.354789
#> 
```
