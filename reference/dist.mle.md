# Maximum likelihood estimation of Spherical Rotational Symmetric Distributions

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

## See also

Other distribution-MLE:
[`bingham-mle`](https://tobiste.github.io/structr/reference/bingham-mle.md),
[`fisher-mle`](https://tobiste.github.io/structr/reference/fisher-mle.md),
[`watson-mle`](https://tobiste.github.io/structr/reference/watson-mle.md)

## Examples

``` r
x <- rkent(100, mu = Line(120, 50), k = 5, b = 1)
kent_mle(x)
#> $G
#> Vector (Vec3) object (n = 3):
#>                x          y           z
#> mean  -0.4050362  0.4805988  0.77779847
#> major -0.8191348 -0.5686510 -0.07519461
#> minor  0.4061574 -0.6675783  0.62399945
#> 
#> $param
#>     kappa      beta       psi 
#> 4.7480627 1.0495075 0.4603117 
#> 
#> $logcon
#> [1] 5.081248
#> 
#> $loglik
#> [1] -129.6656
#> 
vmf_mle(x)
#> $loglik
#> [1] -134.3081
#> 
#> $mu
#> Vector (Vec3) object (n = 1):
#>          x          y          z 
#> -0.4050362  0.4805988  0.7777985 
#> 
#> $kappa
#> [1] 4.452437
#> 
```
