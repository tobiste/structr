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
#> mean  113.8243 57.43571
#> major 340.8044 23.54371
#> minor 241.1014 21.14700
#> 
#> $param
#>      kappa       beta        psi 
#>  4.1174745  0.7651603 -0.4799822 
#> 
#> $logcon
#> [1] 4.574504
#> 
#> $loglik
#> [1] -142.7972
#> 
vmf_mle(x)
#> $loglik
#> [1] -145.5938
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 113.82428  57.43571 
#> 
#> $kappa
#> [1] 3.969948
#> 
```
