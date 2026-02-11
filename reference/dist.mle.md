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
#> mean  110.0929 59.68398
#> major 332.4360 23.37337
#> minor 234.2738 18.18560
#> 
#> $param
#>      kappa       beta        psi 
#>  4.7530889  1.0564822 -0.5983497 
#> 
#> $logcon
#> [1] 5.085824
#> 
#> $loglik
#> [1] -129.5824
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.2927
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 110.09292  59.68398 
#> 
#> $kappa
#> [1] 4.453128
#> 
```
