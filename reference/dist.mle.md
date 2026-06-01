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
#>         azimuth    plunge
#> mean  129.65896 56.664592
#> major  24.31214  9.875323
#> minor 288.19482 31.472646
#> 
#> $param
#>     kappa      beta       psi 
#> 4.1679369 0.8950880 0.2883576 
#> 
#> $logcon
#> [1] 4.623915
#> 
#> $loglik
#> [1] -141.8436
#> 
vmf_mle(x)
#> $loglik
#> [1] -145.767
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 129.65896  56.66459 
#> 
#> $kappa
#> [1] 3.962919
#> 
```
