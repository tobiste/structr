# Bingham test of uniformity of axes on the sphere

Bingham test of uniformity of axes on the sphere

## Usage

``` r
bingham_test(D, n, level = 0.95)
```

## Source

modified after the python package `apsg` (v.2.0.6) by O. Lexa

## Arguments

- D:

  3x3 orientation tensor (see
  [`ortensor()`](https://tobiste.github.io/structr/reference/ortensor.md))

- n:

  integer. Number of vectors

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

## Details

The statistic is \$\$S = 15 n / 2 (\text{tr}(T^2) - 1/3)\$\$ with the
orientation tensor \\T\\ (scaled to trace 1), asymptotically chi-square
with 5 degrees of freedom (MJ eq. 10.7.1-10.7.2), which gives the
p-value.

The modified statistic \\S^\*\\ of MJ eq. 10.7.3, a polynomial
correction meant for \\S\\ near the critical values, is given in
`"extra"`.

\\S\\ equals Bingham's (1974, eq. 5.9) \\X_U^2\\ and \\n I\\ with
Lisle's (1985) intensity \$\$I = 15/2 \sum((\lambda_i - 1/3)^2)\$\$.

Large \\S\\ indicates a preferred orientation (cluster or girdle); a
very small \\S\\ (lower tail) indicates data more regular than random.

## See also

Other tests:
[`fisher_ftest()`](https://tobiste.github.io/structr/reference/fisher_ftest.md),
[`rayleigh_test3()`](https://tobiste.github.io/structr/reference/rayleigh_test3.md),
[`wellner()`](https://tobiste.github.io/structr/reference/wellner.md)

## Examples

``` r
t <- example_lines
bingham_test(ortensor(t), n = nrow(t))
#> Error in bingham_test(ortensor(t), n = nrow(t)): could not find function "bingham_test"
```
