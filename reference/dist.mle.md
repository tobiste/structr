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
#> mean  115.3327 60.83988
#> major 346.5714 19.25591
#> minor 248.8572 21.01942
#> 
#> $param
#>      kappa       beta        psi 
#>  4.2378366  1.0430515 -0.3514472 
#> 
#> $logcon
#> [1] 4.691606
#> 
#> $loglik
#> [1] -140.5019
#> 
vmf_mle(x)
#> $loglik
#> [1] -145.894
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 115.33274  60.83988 
#> 
#> $kappa
#> [1] 3.957774
#> 
```
