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
#>        azimuth     plunge
#> mean  113.8229 50.3988403
#> major 204.5268  0.5823217
#> minor 295.0085 39.5951347
#> 
#> $param
#>     kappa      beta       psi 
#> 5.3599562 1.0318862 0.3438614 
#> 
#> $logcon
#> [1] 5.562864
#> 
#> $loglik
#> [1] -117.8005
#> 
vmf_mle(x)
#> $loglik
#> [1] -121.5514
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 113.82289  50.39884 
#> 
#> $kappa
#> [1] 5.062787
#> 
```
