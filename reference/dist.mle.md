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
#>          azimuth   plunge
#> mean  126.012977 57.55204
#> major   5.086143 18.09542
#> minor 265.938729 25.94369
#> 
#> $param
#>       kappa        beta         psi 
#>  3.87044353  0.67824095 -0.06716341 
#> 
#> $logcon
#> [1] 4.384172
#> 
#> $loglik
#> [1] -148.5554
#> 
vmf_mle(x)
#> $loglik
#> [1] -150.8636
#> 
#> $mu
#> Line object (n = 1):
#>   azimuth    plunge 
#> 126.01298  57.55204 
#> 
#> $kappa
#> [1] 3.76078
#> 
```
