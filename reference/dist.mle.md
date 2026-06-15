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
#> mean  116.2983 57.53567
#> major 341.2430 24.24167
#> minor 241.6921 20.22783
#> 
#> $param
#>      kappa       beta        psi 
#>  4.1027673  0.8213896 -0.4758630 
#> 
#> $logcon
#> [1] 4.568351
#> 
#> $loglik
#> [1] -143.2118
#> 
vmf_mle(x)
#> $loglik
#> [1] -146.5177
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 116.29830  57.53567 
#> 
#> $kappa
#> [1] 3.93258
#> 
```
