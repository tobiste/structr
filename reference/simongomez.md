# Example data set

80 Fault-slip orientations from Simon-Gomez ???

## Usage

``` r
data('simongomez')
```

## Format

An object of class `"Fault`

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
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md)

## Examples

``` r
data("simongomez")
head(simongomez)
#> Fault object (n = 6):
#>      dip_direction dip   azimuth    plunge sense
#> [1,]           310  63  34.01656 11.562272    -1
#> [2,]           253  86 163.20946  2.992685     1
#> [3,]           307  84  35.50710 13.921756    -1
#> [4,]           115  83  26.35703 10.916997     1
#> [5,]           128  85  38.61312  6.973230     1
#> [6,]           238  72 153.06371 15.197478     1
```
