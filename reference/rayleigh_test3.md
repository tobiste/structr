# Rayleigh test of uniformity against a unimodal alternative

Rayleigh test of uniformity against a unimodal alternative

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

## Details

The statistic is \\S = p n \bar{R}^2\\ for \\n\\ unit vectors in \\p\\
dimensions, asymptotically chi-square with \\p\\ degrees of freedom; the
p-value is that of the modified statistic \\S^\* = (1 - 1/(2n)) S + S^2
/ (2n(p + 2))\\, which has error \\O(1/n^2)\\ (MJ eq. 10.4.4-10.4.6, eq.
6.3.3 for the circle; FLE sec. 5.3.1(i)). The test is not consistent
against axial alternatives.

## See also

Other tests:
[`bingham_test()`](https://tobiste.github.io/structr/reference/bingham_test.md),
[`fisher_ftest()`](https://tobiste.github.io/structr/reference/fisher_ftest.md),
[`wellner()`](https://tobiste.github.io/structr/reference/wellner.md)

## Examples

``` r
rayleigh_test3(example_lines)
#> Error in rayleigh_test3(example_lines): could not find function "rayleigh_test3"
```
