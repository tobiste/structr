# Vorticity from Rotated Porphyroclasts

194 measurements of aspect ratio and orientation (wrt. foliation) of
rotated porphyroclasts in a mylonitic, greenschist-grade metavolcanic
rock from the Shebandowan Greenstone Belt (Superior Province, Canada).
Sample HT-17-XZ from Stephan et al. (2025)

## Usage

``` r
data('shebandowan')
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
[`example_planes_df`](https://tobiste.github.io/structr/reference/example_planes_df.md),
[`gray_example`](https://tobiste.github.io/structr/reference/gray_example.md),
[`holst`](https://tobiste.github.io/structr/reference/holst.md),
[`hossack1968`](https://tobiste.github.io/structr/reference/hossack1968.md),
[`ramsay`](https://tobiste.github.io/structr/reference/ramsay.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md)

## Examples

``` r
data("shebandowan")
head(shebandowan)
#>          x        y     r     phi   sense      area
#> 1 2050.462 1896.223 1.542  178.53 dextral 0.4773643
#> 2 1741.010 1575.699 1.092  134.62 dextral 0.8827331
#> 3 3405.425  586.568 1.302   15.35 dextral 0.6845254
#> 4 3354.368 2536.469 1.031   93.22 dextral 1.0000000
#> 5 1595.929 2035.022 1.894  -22.38 dextral 0.1454859
#> 6  912.342  834.246 1.050 -116.57 dextral 0.1397272
```
