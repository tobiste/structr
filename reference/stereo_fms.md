# Plot Focal Mechanism Solution in a Spherical Plot

Shows the principal faults associated with a stress tensor and a given
friction coefficient using the focal mechanism / beach ball style.

## Usage

``` r
stereo_fms(
  sigma,
  friction = 0.6,
  fill = TRUE,
  col = "#BEBEBE80",
  border = "black",
  ...
)
```

## Arguments

- sigma:

  symmetric 3x3 matrix. The (reduced) stress tensor.

- friction:

  numeric. Coefficient of friction (0.6 by default)

- fill:

  logical. Whether to fill the inner part of the small-circle? `FALSE`
  by default.

- col:

  fill color of quadrant

- border:

  Color of the filled small-circle's outline (ignored if `fill=FALSE`)

- ...:

  optional plotting parameters passed to
  [`stereo_greatcircle()`](https://tobiste.github.io/structr/reference/stereo_cones.md)

## Value

`"Fault"` object

## See also

`stereo_fms()` and
[`principal_fault()`](https://tobiste.github.io/structr/reference/principal_fault.md)

## Examples

``` r
f <- angelier1990$TYM
sig <- reduced_stress(f)

stereoplot()
stereo_fms(sig)
```
