# Rayleigh test of uniformity against a unimodal alternative in \\S^3\\

Rayleigh test of uniformity against a unimodal alternative in \\S^3\\

## Usage

``` r
rayleigh_test3(x, level = 0.95)
```

## Source

modified after the python package `apsg` (v.2.0.6) by O. Lexa

## Arguments

- x:

  Object of class `spherical`

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

## Value

object of class `"Test"`

## Details

The statistic is \$\$S = p n \bar{R}^2\$\$ for \\n\\ unit vectors in
\\p\\ dimensions, asymptotically chi-square with \\p\\ degrees of
freedom; the p-value is that of the modified statistic \$\$S^\* = (1 -
\frac{1}{2n}) S + \frac{S^2}{2n(p + 2)}\$\$ which has error \\O(1/n^2)\\
(Mardia and Jupp, 2000; eq. 10.4.4-10.4.6, eq. 6.3.3 for the circle;
Fisher et al., 1987, sec. 5.3.1(i)). The test is not consistent against
axial alternatives.

## References

Mardia, K.V. and Jupp, P.E. (2000) Directional Statistics. John Wiley &
Sons, Chichester.

Fisher, N.I., Lewis, T. and Embleton, B.J.J. (1987) Statistical Analysis
of Spherical Data. Cambridge University Press, Cambridge.

## See also

Other tests:
[`bingham_test()`](https://tobiste.github.io/structr/reference/bingham_test.md),
[`fisher_ftest()`](https://tobiste.github.io/structr/reference/fisher_ftest.md),
[`fisher_goodness_of_fit()`](https://tobiste.github.io/structr/reference/fisher_goodness_of_fit.md),
[`fisher_kent_test()`](https://tobiste.github.io/structr/reference/fisher_kent_test.md),
[`watson_symmetry_test()`](https://tobiste.github.io/structr/reference/watson_symmetry_test.md),
[`wellner()`](https://tobiste.github.io/structr/reference/wellner.md)

## Examples

``` r
rayleigh_test3(example_lines)
#> 
#>  Rayleigh Test
#> 
#> statistic = 1954991, df = 3, p-value < 2.2e-16
#> critical value (95% confidence level) = 7.8147
#> null hypothesis rejected at the 95% confidence level: yes
#> further statistics:
#>            R         Rbar     modified 
#> 7.398632e+03 8.807895e+01 4.551932e+09 
#> source: Mardia & Jupp (2000), eq. 10.4.4-10.4.6 
#> 
```
