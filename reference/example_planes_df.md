# Planes orientations data

200 foliation measurements from a folded metasedimentary sequence in the
Quetico Subprovince (Stephan et al., 2025)

## Usage

``` r
data('example_planes_df')
```

## Format

An object of class `data.frame`

## References

Stephan, T., Phillips, N., Tiitto, H., Perez, A., Nwakanma, M., Creaser,
R., & Hollings, P. (2025). Going with the flow - Changes of vorticity
control gold enrichment in Archean shear zones (Shebandowan Greenstone
Belt, Superior Province, Canada). Journal of Structural Geology, 201,
105542.
[doi:10.1016/j.jsg.2025.105542](https://doi.org/10.1016/j.jsg.2025.105542)

## See also

Other datasets:
[`angelier1990`](https://tobiste.github.io/structr/reference/angelier1990.md),
[`example_lines`](https://tobiste.github.io/structr/reference/example_lines.md),
[`example_lines_df`](https://tobiste.github.io/structr/reference/example_lines_df.md),
[`example_planes`](https://tobiste.github.io/structr/reference/example_planes.md),
[`gray_example`](https://tobiste.github.io/structr/reference/gray_example.md),
[`holst`](https://tobiste.github.io/structr/reference/holst.md),
[`ramsay`](https://tobiste.github.io/structr/reference/ramsay.md),
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md)

## Examples

``` r
data("example_planes_df")
head(example_planes_df)
#> # A tibble: 6 × 4
#>   dipdir   dip quality feature_type
#>    <dbl> <dbl>   <dbl> <chr>       
#> 1    142    52       3 foliation   
#> 2    135    43       3 foliation   
#> 3    148    42       3 foliation   
#> 4    150    46       3 foliation   
#> 5    139    51       3 foliation   
#> 6    158    51       3 foliation   
```
