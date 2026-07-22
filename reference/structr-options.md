# Global projection options for structr

Plotting defaults for stereographic/equal-area projections. For example,
`structr.earea` and `structr.upper.hem` control the default hemisphere
projection used across plotting and labeling functions
([`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md),
[`plot.Vec3()`](https://tobiste.github.io/structr/reference/plot-spherical.md),
[`text.spherical()`](https://tobiste.github.io/structr/reference/text.spherical.md),
...).

Set them with [`options()`](https://rdrr.io/r/base/options.html) to
change the default for an entire session; any function argument passed
explicitly overrides the global default for that call.

- `structr.earea`:

  Logical. Equal-area (`TRUE`, default) vs. equal-angle (`FALSE`)
  projection.

- `structr.upper.hem`:

  Logical. Upper (`TRUE`) vs. lower hemisphere (`FALSE`, default).

- `structr.guides`:

  Logical. logical. Whether guides should be added to the plot (`TRUE`
  by default).

- `structr.d`:

  integer. Angle distance between guides. Default: `10`

- `structr.col`:

  Color of guide lines. `"gray90"` by default.

- `structr.lwd`:

  Width of guide lines. `0.5` by default.

- `structr.lty`:

  Type of guide lines. `1` by default.

- `structr.border.col`:

  Color of primitive circle (frame), center-cross and ticks of the
  stereo plot. `"black"` by default.

- `structr.centercross`:

  Logical. Whether a center cross should be added (`TRUE` by default).

- `structr.ticks`:

  Integer. Angle between ticks. if `NULL` (the default), no ticks are
  drawn.

- `structr.origin.text`:

  character. Text at origin of plot. `"N"` by default.

- `structr.labels`:

  this can either be a logical value specifying whether (numerical)
  annotations are to be made next to the tick marks, or a character or
  expression

- `structr.ladj`:

  adjustment for all labels away from origin of projection circle. This
  essentially an amount that is added to `radius` and the length of the
  ticks.

- `structr.radius`:

  Numeric. Radius of the projection circle. `1` by default.

## See also

[`stereoplot()`](https://tobiste.github.io/structr/reference/stereoplot.md)

## Examples

``` r
old <- options(structr.upper.hem = TRUE)
# ... plotting calls now default to upper hemisphere
options(old)  # restore
```
