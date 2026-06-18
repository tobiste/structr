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
#> mean  125.8265 57.775859
#> major 217.8369  1.266775
#> minor 308.6346 32.193094
#> 
#> $param
#>     kappa      beta       psi 
#> 4.0874759 0.5846227 0.5897405 
#> 
#> $logcon
#> [1] 4.538854
#> 
#> $loglik
#> [1] -143.313
#> 
vmf_mle(x)
#> $loglik
#> [1] -144.8093
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 125.82647  57.77586 
#> 
#> $kappa
#> [1] 4.001921
#> 
```
