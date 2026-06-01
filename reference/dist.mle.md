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
#> mean  120.9238 52.52866
#> major  27.8533  2.35119
#> minor 296.0562 37.37143
#> 
#> $param
#>     kappa      beta       psi 
#> 4.4488711 0.9865786 0.3763238 
#> 
#> $logcon
#> [1] 4.845038
#> 
#> $loglik
#> [1] -135.7829
#> 
vmf_mle(x)
#> $loglik
#> [1] -140.2288
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 120.92380  52.52866 
#> 
#> $kappa
#> [1] 4.193136
#> 
```
