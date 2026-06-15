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
#> mean  126.83740 53.255748
#> major  31.53223  3.948812
#> minor 298.60858 36.460423
#> 
#> $param
#>     kappa      beta       psi 
#> 4.3063891 0.4775189 0.4252541 
#> 
#> $logcon
#> [1] 4.698845
#> 
#> $loglik
#> [1] -138.2025
#> 
vmf_mle(x)
#> $loglik
#> [1] -139.0097
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 126.83740  53.25575 
#> 
#> $kappa
#> [1] 4.245379
#> 
```
