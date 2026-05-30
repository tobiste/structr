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
#> mean  117.46504 54.839239
#> major  23.75737  2.608052
#> minor 291.92726 35.034669
#> 
#> $param
#>     kappa      beta       psi 
#> 4.9775378 1.3916597 0.3227297 
#> 
#> $logcon
#> [1] 5.297296
#> 
#> $loglik
#> [1] -126.0612
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.1353
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 117.46504  54.83924 
#> 
#> $kappa
#> [1] 4.46022
#> 
```
