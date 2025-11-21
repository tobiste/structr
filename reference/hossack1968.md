# 3D Strain Data II

Example data from J.R. Hossack (1968) containing 3D strain data from 81
localities from the Bygdin Conglomerate by the Upper Jotun Nappe. The
data gives the mean axial ratios (X, Y, Z) of the deformed pebbles.
Different to Hossack, here X is the longest axes, and Z the shortest.

## Usage

``` r
data('hossack1968')
```

## Format

An object of class `matrix`

## See also

Other datasets:
[`angelier1990`](https://tobiste.github.io/structr/reference/angelier1990.md),
[`example_lines`](https://tobiste.github.io/structr/reference/example_lines.md),
[`example_lines_df`](https://tobiste.github.io/structr/reference/example_lines_df.md),
[`example_planes`](https://tobiste.github.io/structr/reference/example_planes.md),
[`example_planes_df`](https://tobiste.github.io/structr/reference/example_planes_df.md),
[`gray_example`](https://tobiste.github.io/structr/reference/gray_example.md),
[`holst`](https://tobiste.github.io/structr/reference/holst.md),
[`ramsay`](https://tobiste.github.io/structr/reference/ramsay.md),
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md)

## Examples

``` r
data("hossack1968")
head(hossack1968)
#>   Z    Y    X
#> 1 1 11.1 18.3
#> 2 1  4.2 32.0
#> 3 1  2.3 30.0
#> 4 1  4.0 16.6
#> 5 1  7.5 54.3
#> 6 1  7.4 30.2
```
