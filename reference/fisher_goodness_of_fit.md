# Goodness-of-fit tests of the Fisher distribution

Goodness-of-fit tests of the Fisher distribution

## Usage

``` r
fisher_goodness_of_fit(x, level = 0.95)
```

## Source

modified after the python package `apsg` (v.2.0.6) by O. Lexa

## Arguments

- x:

  object of class `"Spherical"`

- level:

  numeric. Confidence level (\\1-\alpha\\). Default is `0.95`

## Value

object of class `"Test"`

## Details

The data are rotated to their mean direction. Under the Fisher model the
colatitudes \\\theta'\\ give \\1 - \cos(\theta')\\ approximately
exponential with rate \\k = (n - 1) / \sum(1 - \cos(\theta'))\\ and the
longitudes are uniform.

The colatitude test \$\$M_E = (D_n - 0.2/n)(\sqrt(n) + 0.26 +
0.5/\sqrt(n))\$\$ uses the Kolmogorov statistic \\D_n\\ against the
exponential distribution, while the longitude test \$\$M_U =
V_n(\sqrt(n) - 0.467 + 1.623/\sqrt(n))\$\$ uses the Kuiper statistic
\\V_n\\ against the uniform distribution (FLE eq. 5.15-5.16).

Critical values are known for the levels 0.90, 0.95 and 0.99.

## Note

The constants are FLE's, which are the original ones of Stephens (1974)
and Fisher & Best (1984). Mardia & Jupp (2000, eq. 12.3.1, 12.3.4)
misquote them.

## Examples

``` r
set.seed(20250411)
x <- rvmf(100)
fisher_goodness_of_fit(x)
#> Error in fisher_goodness_of_fit(x): could not find function "fisher_goodness_of_fit"
```
