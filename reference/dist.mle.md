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
#> mean  121.7504 55.979595
#> major 223.3169  7.707913
#> minor 318.3401 32.899967
#> 
#> $param
#>    kappa     beta      psi 
#> 4.548943 1.019742 0.716004 
#> 
#> $logcon
#> [1] 4.924794
#> 
#> $loglik
#> [1] -133.7402
#> 
vmf_mle(x)
#> $loglik
#> [1] -138.3786
#> 
#> $mu
#> Line object (n = 1):
#>  azimuth   plunge 
#> 121.7504  55.9796 
#> 
#> $kappa
#> [1] 4.272651
#> 
```
