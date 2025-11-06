# Intersection lineation data

84 intersection lineations from a folded metasedimentary sequence in the
Quetico Subprovince (Stephan et al., 2025)

## Usage

``` r
data('example_lines_df')
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
[`example_planes`](https://tobiste.github.io/structr/reference/example_planes.md),
[`example_planes_df`](https://tobiste.github.io/structr/reference/example_planes_df.md),
[`gray_example`](https://tobiste.github.io/structr/reference/gray_example.md),
[`holst`](https://tobiste.github.io/structr/reference/holst.md),
[`ramsay`](https://tobiste.github.io/structr/reference/ramsay.md),
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md)

## Examples

``` r
data("example_lines_df")
head(example_lines_df)
#> # A tibble: 6 × 4
#>   trend plunge quality feature_type
#>   <dbl>  <dbl>   <dbl> <chr>       
#> 1    54     13       3 stretching  
#> 2    61     15       3 stretching  
#> 3    74     14      NA stretching  
#> 4    80     19      NA stretching  
#> 5    63     17      NA stretching  
#> 6    76     10      NA stretching  
```
