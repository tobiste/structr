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
#> mean  110.1879 52.027091
#> major  19.7347  0.353697
#> minor 289.4586 37.970658
#> 
#> $param
#>     kappa      beta       psi 
#> 4.6578446 0.7975418 0.2720851 
#> 
#> $logcon
#> [1] 4.989649
#> 
#> $loglik
#> [1] -130.9964
#> 
vmf_mle(x)
#> $loglik
#> [1] -133.5775
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 110.18786  52.02709 
#> 
#> $kappa
#> [1] 4.485434
#> 
```
