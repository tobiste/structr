# 3D Strain Data

Example data from Holst and Fossen (1987) containing 3D strain data from
17 localities in the West Norwegian Caledonides. The data gives the
difference of the \\\epsilon_1\\ and \\\epsilon_2\\ (`e1e2`) and
\\\epsilon_2\\ and \\\epsilon_3\\ (`e2e3`)

## Usage

``` r
data('holst')
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
[`ramsay`](https://tobiste.github.io/structr/reference/ramsay.md),
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md)

## Examples

``` r
data("holst")
head(holst)
#>          R_XY      R_YZ
#> [1,] 1.030455 30.876643
#> [2,] 1.072508 27.112639
#> [3,] 1.040811 25.028120
#> [4,] 1.072508 27.938342
#> [5,] 1.150274 20.085537
#> [6,] 2.117000  9.025013
```
