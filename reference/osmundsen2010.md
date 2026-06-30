# Fault-Slip Data from Northern Norway

198 measurements fault-slip measurements from a quarry in Røsand on the
Lofoten-Vesterålen archipelago in Northern Norway

## Usage

``` r
data('osmundsen2010')
```

## Format

An object of class `Pair`

## References

Osmundsen, P.T., Redfield, T.F., Hendriks, B.H.W., Bergh, S.G., Hansen,
J.-A., Henderson, I.H.C., Dehls, J., Lauknes, T.R., Larsen, Y., Anda,
E., Davidsen, B., 2010. Fault-controlled Alpine topography in Norway.
Journal of the Geological Society 167, 83-98.

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
[`shebandowan`](https://tobiste.github.io/structr/reference/shebandowan.md),
[`simongomez`](https://tobiste.github.io/structr/reference/simongomez.md),
[`strabo_prj`](https://tobiste.github.io/structr/reference/strabo_prj.md)

## Examples

``` r
data("osmundsen2010")
head(osmundsen2010)
#> Pair object (n = 6):
#>   dip_direction dip    azimuth plunge
#> 1           312  52 312.000000     52
#> 2           310  64 310.000000     64
#> 3           314  40 340.096982     37
#> 4           310  48 310.000000     48
#> 5           348  47 332.938540     46
#> 6           342  48   3.188608     46
```
