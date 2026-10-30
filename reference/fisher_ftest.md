# Test of Mean Orientations

Test against the null-hypothesis that the samples are drawn from the
same Fisher population.

## Usage

``` r
fisher_ftest(x, y, level = 0.95, na.rm = TRUE)
```

## Arguments

- x, y:

  objects of class `"Vec3"`, `"Line"`, `"Ray"`, or `"Plane"`, where the
  rows are the observations and the columns are the coordinates.

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

- na.rm:

  logical. Whether `NA` values should be removed before the computation
  proceeds.

## Value

object of class `"Test"`

## See also

Other tests:
[`bingham_test()`](https://tobiste.github.io/structr/reference/bingham_test.md),
[`fisher_goodness_of_fit()`](https://tobiste.github.io/structr/reference/fisher_goodness_of_fit.md),
[`fisher_kent_test()`](https://tobiste.github.io/structr/reference/fisher_kent_test.md),
[`rayleigh_test3()`](https://tobiste.github.io/structr/reference/rayleigh_test3.md),
[`watson_symmetry_test()`](https://tobiste.github.io/structr/reference/watson_symmetry_test.md),
[`wellner()`](https://tobiste.github.io/structr/reference/wellner.md)

## Examples

``` r
set.seed(20250411)
x <- rvmf(100, mu = Line(120, 50), k = 20)
y <- rvmf(100, mu = Line(180, 45), k = 20)

stereoplot()
stereo_point(x, col = 1)
stereo_point(y, col = 2)


fisher_ftest(x, y)
#> 
#>  Fisher F-test for a common mean direction
#> 
#> statistic = 237.14, df = 2, df = 396, p-value < 2.2e-16
#> critical value (95% confidence level) = 3.0185
#> null hypothesis rejected at the 95% confidence level: yes
#> further statistics:
#>        Rx        Ry         R 
#>  95.38633  95.35390 179.65004 
#> source: Fisher, Lewis & Embleton (1987), sec. 7.2.2 
#> 
```
