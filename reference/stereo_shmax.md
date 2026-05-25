# Horizontal directions

Horizontal directions

## Usage

``` r
stereo_shmax(
  azi,
  ...,
  type = c("arrows", "line"),
  shmin = FALSE,
  BALL.radius = 1,
  arrow.offset = 0.02,
  arrow.length = 0.1,
  arrow.head = arrow.length
)
```

## Arguments

- azi:

  numeric. Angle of maximum horizontal stress in degrees

- ...:

  ... arguments passed to
  [`graphics::segments()`](https://rdrr.io/r/graphics/segments.html) or
  [`graphics::arrows()`](https://rdrr.io/r/graphics/arrows.html)

- type:

  character. Either `'line'` or `'arrows'` for a straight line or arrows
  at the perimeter of the plot.

- shmin:

  logical. Whether the minimum horizontal stress should be indicated
  too?

- BALL.radius:

  numeric. Radius of the stereo plot.

- arrow.offset:

  numeric. Offset of the arrows from the perimeter.

- arrow.length:

  numeric. Length of the arrows.

- arrow.head:

  numeric. Length of the arrow head.

## Examples

``` r
stereoplot()
stereo_shmax(30, shmin = TRUE)
```
