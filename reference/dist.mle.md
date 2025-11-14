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
#> mean  117.3363 49.52229
#> major  18.7143  7.29074
#> minor 282.6516 39.54075
#> 
#> $param
#>     kappa      beta       psi 
#> 5.0015031 1.0686473 0.1778848 
#> 
#> $logcon
#> [1] 5.281088
#> 
#> $loglik
#> [1] -124.6798
#> 
vmf_mle(x)
#> $loglik
#> [1] -129.1592
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 117.33631  49.52229 
#> 
#> $kappa
#> [1] 4.689888
#> 
```
