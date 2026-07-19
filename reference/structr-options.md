# Global projection options for yourpkg

`structr.earea` and `structr.upper.hem` control the default hemisphere
projection used across plotting and labeling functions
([`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`plot.Vec3()`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md),
...). Set them with [`options()`](https://rdrr.io/r/base/options.html)
to change the default for an entire session; any function argument
passed explicitly overrides the global default for that call.

- `structr.earea`:

  Logical. Equal-area (`TRUE`, default) vs. equal-angle (`FALSE`)
  projection.

- `structr.upper.hem`:

  Logical. Upper (`TRUE`) vs. lower hemisphere (`FALSE`, default).

- `structr.guides`:

  Logical. logical. Whether guides should be added to the plot (`TRUE`
  by default)

- `structr.radius`:

  numeric. Radius of the projection circle. `1` by default.

## Examples

``` r
old <- options(structr.upper.hem = TRUE)
# ... plotting calls now default to upper hemisphere
options(old)  # restore
```
